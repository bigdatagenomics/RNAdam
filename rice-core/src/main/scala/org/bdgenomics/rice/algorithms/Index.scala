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
package org.bdgenomics.rice.algorithms

import org.apache.spark.Logging
import org.apache.spark.rdd.RDD
import org.apache.spark.SparkContext._
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.formats.avro.NucleotideContigFragment
import org.bdgenomics.rice.Timers._
import org.bdgenomics.rice.models.{ Exon, Transcript }
import org.bdgenomics.adam.util.{ TwoBitFile, ReferenceFile }
import scala.util.hashing.MurmurHash3

object Index extends Serializable with Logging {

  /**
   * Computes an index, given a reference sequence, a set of transcripts, and a k-mer
   * length. An index provides bidirectional mappings between k-mers and equivalence classes.
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
   *
   * @param referenceFile A ReferenceFile representing the chromosome
   * @param transcripts An RDD containing transcripts.
   * @param kmerLength The length _k_ of the k-mers for us to index.
   * @return Returns a tuple containing two RDDs: the first RDD maps k-mers to equivalence class
   *         IDs, the second maps equivalence class IDs to an iterable of member k-mers.
   */
  def apply(referenceFile: ReferenceFile,
            transcripts: RDD[Transcript],
            kmerLength: Int): (RDD[(String, Long)], RDD[(Long, Iterable[String])]) = {

    findEquivalenceClasses(transcripts, kmerLength, referenceFile)
  }

  /**
   * Given an RDD of transcripts, this method finds the k-mer equivalence classes. K-mer
   * equivalence classes represent k-mers that show up with equal abundance across several
   * transcripts---in our method, we simplify by just looking for exons that occur across
   * several transcripts. This method returns two RDDs: one maps equivalence class IDs to
   * a list of transcripts, and the other maps exon IDs to equivalence class IDs.
   *
   * @param transcripts An RDD of transcripts.
   * @param kmerLength The length of kmers to calssify by
   * @param referenceFile The ReferenceFile representing the chromosome
   * @return Returns two RDDs: one maps equivalence class IDs to a list of transcripts,
   *         and the other maps exon IDs to equivalence class IDs.
   */
  private[algorithms] def findEquivalenceClasses(transcripts: RDD[Transcript],
                                                 kmerLength: Int,
                                                 referenceFile: ReferenceFile): (RDD[(String, Long)], RDD[(Long, Iterable[String])]) = {
    // Broadcast variable representing the reference file:
    val refFile = Broadcast.time {
      transcripts.context.broadcast(referenceFile)
    }

    // RDD of (list of kmers in eq class, equivalence class ID)
    val kmersToClasses = GenerateClasses.time {
      val kmersAndTranscript = KmersAndTranscript.time {
        transcripts.flatMap(t => {
          val sequence = Extract.time {
            refFile.value.extract(t.region)
          }
          val kmers = SplitKmers.time {
            sequence.sliding(kmerLength).toIterable
          }
          kmers.map(k => ((t.id, k), 1))
        })
      }

      val kmersByCount = CollectingKmersByCount.time { kmersAndTranscript.reduceByKey(_ + _) }
      val sortByTranscript = SortByTranscript.time { kmersByCount.map(v => ((v._1._1, v._2), v._1._2)) }
      val kmersByTranscript = CollectingKmersByTranscript.time { sortByTranscript.groupByKey() }
      val eqClasses = DistillingEqClasses.time { kmersByTranscript.map(v => v._2) }
      val numberedEqClasses = NumberingEqClasses.time { eqClasses.zipWithUniqueId() }
      numberedEqClasses
    }

    GenerateIndices.time {
      // RDD of (kmer, equivalence class ID)
      val kmerToClassID = kmersToClasses.flatMap(v => {
        v._1.map((_, v._2))
      })

      // RDD of (equivalence class ID, list of kmers)
      val idsToKmers = kmersToClasses.map(v => (v._2, v._1))

      (kmerToClassID, idsToKmers)
    }
  }
}
