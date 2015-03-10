/**
 * Licensed to Big Data Genomics (BDG) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The BDG licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.bdgenomics.RNAdam.algorithms.quantification

import org.apache.spark.Logging
import org.apache.spark.rdd.RDD
import org.apache.spark.SparkContext._
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.formats.avro.AlignmentRecord
import org.bdgenomics.adam.models.Transcript

object Quantify extends Serializable with Logging {

  /**
   *
   * This code is based on the implementation of Sailfish, which is described in:
   *
   * Patro, Rob, Stephen M. Mount, and Carl Kingsford. "Sailfish: alignment-free isoform
   * quantification from RNA-seq reads using lightweight algorithms." arXiv preprint arXiv:1308.3700 (2013).
   *
   * and:
   *
   * Patro, Rob, Stephen M. Mount, and Carl Kingsford. "Sailfish enables alignment-free
   * isoform quantification from RNA-seq reads using lightweight algorithms." Nature biotechnology 32.5 (2014): 462-464.
   */
  def apply(reads: RDD[AlignmentRecord],
            kmerToEquivalenceClass: RDD[(String, Long)],
            equivalenceClassToTranscript: RDD[(Long, Iterable[String])],
            transcripts: RDD[Transcript],
            kmerLength: Int,
            maxIterations: Int,
            calibrateKmerBias: Boolean = true,
            calibrateLengthBias: Boolean = true): RDD[(Transcript, Double)] = {

    // cache transcripts, then compute transcript lengths
    transcripts.cache()
    val tLen = extractTranscriptLengths(transcripts)

    // cut reads into kmers and then calibrate if desired
    val readKmers = reads.adamCountKmers(kmerLength)
    val calibratedKmers = if (calibrateKmerBias) {
      Tare.calibrateKmers(readKmers)
    } else {
      readKmers
    }

    // map kmer counts into equivalence classes
    val equivalenceClassCounts = mapKmersToClasses(calibratedKmers, kmerToEquivalenceClass)

    // Cache the RDD equivalenceClassCounts so that it is not computed twice.
    equivalenceClassCounts.cache()

    // The total of all the equivalence class counts.
    // Should be equal to the number of kmers.
    val numKmers: Long = equivalenceClassCounts.map(kv => kv._2).reduce(_ + _)

    // The relative number of kmers in each equivalence class.
    // This is needed by the maximization step of the EM algorithm below
    val relNumKmersInEC = equivalenceClassCounts.map((ec: (Long, Long)) => {
      (ec._1, ec._2.toDouble / numKmers)
    }).collectAsMap

    // we initialize the alphas by splitting all counts equally across transcripts
    var alpha = initializeEM(equivalenceClassCounts,
      equivalenceClassToTranscript)

    // we initialize the µ-hat by running a first step of the M algorithm
    var muHat = m(alpha, tLen, kmerLength, relNumKmersInEC)

    // run iterations of the em algorithm
    (0 until maxIterations).foreach(i => {
      log.info("On iteration " + i + " of EM algorithm.")

      alpha = e(muHat)
      muHat = m(alpha, tLen, kmerLength, relNumKmersInEC)
    })

    // perform calibration to correct for transcript length bias, if desired
    if (calibrateLengthBias) {
      muHat = Tare.calibrateTxLenBias(muHat, tLen)
    }

    // join transcripts up and return
    joinTranscripts(transcripts, muHat)
  }

  /**
   * Extracts transcripts from an RDD of transcripts. Each transcript has an iterable of
   * exons, each of which contains the region covered by the exon. We can derive length
   * from +-reducing the length of these exons.
   *
   * @param rdd RDD of transcripts.
   * @return Returns a map from transcript ID to the transcript length.
   */
  private[quantification] def extractTranscriptLengths(rdd: RDD[Transcript]): scala.collection.Map[String, Long] = {
    rdd.map(t => {
      (t.id, t.exons.map(_.region.width - 1).reduce(_ + _))
    }).collectAsMap
  }

  /**
   * This function takes in a set of counted kmers and a mapping between kmers and
   * equivalence classes and returns the total count of kmers per each equivalence
   * class.
   *
   * @param kmerCounts An RDD of tuples containing (the k-mer string, the number seen).
   * @param kmerToEquivalenceClass An RDD containing tuples which map k-mer strings to
   *                               equivalence class IDs.
   * @return Returns an RDD containing tuples with (equivalence class ID, count).
   */
  private[quantification] def mapKmersToClasses(kmerCounts: RDD[(String, Long)],
                                                kmerToEquivalenceClass: RDD[(String, Long)]): RDD[(Long, Long)] = {
    kmerToEquivalenceClass.join(kmerCounts)
      .map((x: (String, (Long, Long))) => x._2)
      .reduceByKey((c0: Long, c1: Long) => c0 + c1)
  }

  /**
   * Initializes the EM loop by splitting all equivalence class counts equally across
   * all transcripts. Takes in the total coverage count for all equivalence classes and
   * the equivalence class to transcript mapping, and returns an RDD containing a tuple
   * indicating the (transcript ID, normalized coverage, and equivalence classes that
   * map to the transcript).
   *
   * @param equivalenceClassCounts An RDD containing tuples of (each equivalence class ID,
   *                                                            the equivalence class coverage).
   * @param equivalenceClassToTranscript An RDD of tuples mapping equivalence class IDs to
   *                                     transcript IDs.
   * @return Returns an RDD containing tuples of (transcript ID,
   *                                              normalized coverage,
   *                                              iterable of equivalence class IDs).
   */
  private[quantification] def initializeEM(equivalenceClassCounts: RDD[(Long, Long)],
                                           equivalenceClassToTranscript: RDD[(Long, Iterable[String])]): RDD[(Long, Iterable[(String, Double)])] = {
    equivalenceClassCounts.join(equivalenceClassToTranscript).map((x: (Long, (Long, Iterable[String]))) => {
      val normCoverage: Double = x._2._1.toDouble / x._2._2.size()
      val iter2: Iterable[(String, Double)] = x._2._2.map((y: String) => {
        (y, normCoverage)
      })
      (x._1, iter2)
    })
  }

  /**
   * The expectation stage assigns a fraction of the count total per equivalence class to
   * each transcript that the equivalence class belongs to.
   *
   * Per equivalence class s_j and transcript t_i, the update is:
   *
   * α(j,i) = \frac{µhat_i T(s_j)}{\sum_{t \supseteq s_j} µhat_t}
   *
   * @param transcriptWeights An RDD of tuples where each tuple contains a transcript ID,
   *                          the normalized coverage for that transcript, and an iterable
   *                          over the IDs of the equivalency classes that make up this transcript.
   * @return Returns an RDD of tuples which map equivalence class IDs to an iterable of tuples which
   *         map transcript IDs to alpha assignments.
   */
  private[quantification] def e(transcriptWeights: RDD[(String, Double, Iterable[Long])]): RDD[(Long, Iterable[(String, Double)])] = {
    transcriptWeights.flatMap((x: (String, Double, Iterable[Long])) => {
      x._3.map((y: Long) => (y, (x._1, x._2)))
    }).groupByKey()
      .map((c: (Long, Iterable[(String, Double)])) => {
        val tot: Double = c._2.map((tw: (String, Double)) => {
          tw._2
        }).reduce((w0: Double, w1: Double) => {
          w0 + w1
        })
        (c._1, c._2.map((p: (String, Double)) => (p._1, p._2 / tot)))
      })
  }

  /**
   * The maximization step computes the relative abundance of each transcript, given the alpha
   * assignments from each equivalence class to each transcript.
   *
   * Per transcript, we first perform an update:
   *
   * µ_i = \frac{\sum_{s_j \subseteq t_i} α(j,i)}{ * k_j}{lhat_i}
   *
   * k_j is the relative number of k-mers in equivalence class s_j.
   * lhat_i is the adjusted length of transcript i, which equals lhat_i = l_i - k + 1.
   *
   * We then normalize all µ_i by:
   *
   * µhat_i = \frac{µ_i}{\sum_{t_j \in T} µ_j}
   *
   * @param equivalenceClassAssignments An RDD of tuples which map equivalence class IDs to an iterable of
   *                                    tuples which map transcript IDs to alpha assignments.
   * @param tLen A map assigning transcript IDs to the transcript length.
   * @param kmerLength The length of the k-mers used for quantification. Used to adjust the transcript length.
   * @param relNumTxInEC The relative number of transcripts in each equivalence class.
   * @return Returns an RDD containing tuples of (transcript ID,
   *                                              normalized coverage,
   *                                              iterable of equivalence class IDs).
   */
  private[quantification] def m(equivalenceClassAssignments: RDD[(Long, Iterable[(String, Double)])],
                                tLen: scala.collection.Map[String, Long],
                                kmerLength: Int,
                                relNumKmersInEC: scala.collection.Map[Long, Double]): RDD[(String, Double, Iterable[Long])] = {

    // broadcast the transcript length map
    val tLenBcast = equivalenceClassAssignments.context.broadcast(tLen)

    val mus: RDD[(String, Double, Iterable[Long])] = equivalenceClassAssignments.flatMap((eca: (Long, Iterable[(String, Double)])) => {
      eca._2.map((ta: (String, Double)) => {
        (ta._1, (ta._2 * relNumKmersInEC(eca._1), eca._1))
      })
    }).groupByKey()
      .map((raw: (String, Iterable[(Double, Long)])) => {
        val sum_alpha: Double = raw._2.reduce((ae0: (Double, Long), ae1: (Double, Long)) => {
          (ae0._1 + ae1._1, ae0._2)
        })._1
        val adj_len: Long = (tLen(raw._1) - kmerLength) + 1
        val mu: Double = sum_alpha / adj_len
        (raw._1, mu, raw._2.map((ae: (Double, Long)) => ae._2))
      })

    // Cache the mus RDD so that it is not computed twice.
    mus.cache()

    // This normalizes all the µ_i
    // µhat_i = \frac{µ_i}{\sum_{t_j \in T} µ_j}
    val total_mu: Double = mus.map((tme: (String, Double, Iterable[Long])) => {
      tme._2
    }).reduce((m0: Double, m1: Double) => {
      m0 + m1
    })

    // Returns the normalized result.
    mus.map((tme: (String, Double, Iterable[Long])) => {
      (tme._1, tme._2 / total_mu, tme._3)
    })
  }

  /**
   * This method joins the final normalized transcript expression weights against an RDD containing
   * the full descriptors for transcripts.
   *
   * @param transcripts An RDD containing the full, detailed transcript descriptors.
   * @param transcriptWeights An RDD containing tuples of (transcript ID, normalized coverage, and
   *                          an iterable of equivalence class IDs).
   * @return Returns an RDD containing tuples of the full transcript descriptor and the normalized coverage.
   */
  private[quantification] def joinTranscripts(transcripts: RDD[Transcript],
                                              transcriptWeights: RDD[(String, Double, Iterable[Long])]): RDD[(Transcript, Double)] = {

    // RDD of ( ID, coverage ):
    val coverage = transcriptWeights.map(t => (t._1, t._2))

    // Goes from RDD of (Transcript) -> RDD of (ID, Transcript) 
    //-> RDD of (ID, (Transcript, Coverage)) -> RDD of (Transcript, Coverage)
    transcripts.keyBy(t => t.id).join(coverage).map(t => t._2)
  }
}
