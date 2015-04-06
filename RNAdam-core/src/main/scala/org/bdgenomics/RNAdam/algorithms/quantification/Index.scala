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
import org.bdgenomics.formats.avro.NucleotideContigFragment
import org.bdgenomics.adam.models.{ Exon, Transcript }
import org.bdgenomics.adam.util.{ TwoBitFile, ReferenceFile }

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
   * Given a transcript, length and 2BitFile, this method extracts kmers of the specified 
   * length from the transcript
   *
   *@param transcript A Transcript object
   *@param kmerLength The length of kmers to extract
   *@param referenceFile The ReferenceFile representing the chromosome
   *@return Returns a list of kmers contained in Transcript
   *
   */
  private[quantification] def extractKmers(transcript: Transcript,
                                           kmerLength: Int, 
                                           referenceFile: ReferenceFile): Iterable[String] = {
    referenceFile.extract(transcript.region).sliding(kmerLength).toIterable
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
  private[quantification] def findEquivalenceClasses(transcripts: RDD[Transcript],
                                                     kmerLength: Int,
                                                     referenceFile: ReferenceFile): (RDD[(String, Long)], RDD[(Long, Iterable[String])]) = {
    // Broadcast variable representing the reference file:
    val refFile = transcripts.context.broadcast(referenceFile)

    // RDD of (list of kmers in eq class, equivalence class ID)
    val kmersToClasses = transcripts.flatMap(t => {for(k <- extractKmers(t, kmerLength, refFile.value)) 
                                                  yield ((t.id, k), 1)} )                       // ((t.id, kmer), 1)
                                    .foldByKey(0)(_+_)                                          // ((t,id, kmer), abundance)
                                    .map(v => ((v._1._1, v._2), v._1._2))                       // ((t.id, abundance), kmer)
                                    .groupByKey()                                               // ((t.id, abundance), all kmers with same pattern) = (eq class, all kmers in class)
                                    .map(v => v._2)                                             // (kmers in eq class)
                                    .zipWithIndex()                                             // (kmers in eq class, class ID)

    // RDD of (kmer, equivalence class ID)
    val kmerToClassID = kmersToClasses.flatMap(v => {for(k <- v._1) yield (k, v._2)})           // (kmer, class ID)

    // RDD of (equivalence class ID, list of kmers)
    val idsToKmers = kmersToClasses.map(v => (v._2, v._1))                                      // (class ID, list of kmers)

    (kmerToClassID, idsToKmers)
  }

}