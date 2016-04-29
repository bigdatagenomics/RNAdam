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

import org.apache.spark.rdd.RDD
import org.bdgenomics.formats.avro.{ Contig, NucleotideContigFragment, AlignmentRecord }
import org.bdgenomics.adam.models.{ Exon, ReferenceRegion, Transcript, CDS, UTR }
import org.bdgenomics.rice.models.{ KmerIndex, IndexMap }
import org.bdgenomics.rice.algorithms.alignment.{ AlignmentModel, SimpleAlignmentModel }
import org.bdgenomics.rice.utils.riceFunSuite
import net.fnothaft.ananas.models._
import net.fnothaft.ananas.avro.{ Kmer, Backing }
import net.fnothaft.ananas.debruijn.ColoredDeBruijnGraph

object TestAlignmentModel extends AlignmentModel {

  def processRead(iter: Iterator[CanonicalKmer],
                  kmerIndex: KmerIndex): Map[String, Double] = {

    iter.toList
      .flatMap(c => kmerIndex.getTranscripts(c))
      .groupBy(_._1) // Map[ TranscriptID -> List[(TranscriptId, Occurrences)] ] 
      .map(v => v._2.reduce((a, b) => (a._1, a._2 + b._2)))
      .map(v => (v._1, v._2.toDouble))
  }
}

class QuantifySuite extends riceFunSuite {

  // Blatant copy of function from ananas/models/ContigFragment to get around protected status
  def buildFromNCF(fragment: NucleotideContigFragment): ContigFragment = {
    // is this the last fragment in a contig?
    val isLast = fragment.getFragmentNumber == (fragment.getNumberOfFragmentsInContig - 1)

    val sequence = IntMer.fromSequence(fragment.getFragmentSequence)
      .map(_.asInstanceOf[CanonicalKmer])

    new ContigFragment(fragment.getContig.getContigName,
      sequence,
      isLast,
      Option(fragment.getFragmentStartPosition).fold(0)(_.toInt))
  }

  /**
   * Creates a fake contig fragment from the specified sequence.
   *
   * @param sequence A list of Strings which describe the base sequences of the contig fragments we are 
   *                  emulating
   * @param name What to call the resulting contig fragment
   * @return A ContigFragment object that represents the sequence
   */
  def createContigFragment(sequence: String, name: String): ContigFragment = {
    val ncf = NucleotideContigFragment.newBuilder()
      .setContig(Contig.newBuilder()
        .setContigName(name)
        .build())
      .setFragmentNumber(0)
      .setNumberOfFragmentsInContig(1)
      .setFragmentSequence(sequence)
      .build()

    buildFromNCF(ncf)
  }

  /**
   * Computes the correct results for a run of Index upon a set of fragments described by the 
   * sequences and transcripts given.
   *
   * @param sequences A list of Strings which describe the base sequences of the contig fragments we are 
   *                  emulating
   * @param transcriptNames A list of the names of the transcripts corresponding to the sequences in the 
   *                    previous parameter.
   * @return A manually computed index result
   */
  def expectedResults(sequences: Array[String], transcriptNames: Array[String], kmerLength: Int = 16): Map[(Long, Boolean, String), Map[String, Long]] = {
    val combined = sequences.zip(transcriptNames) // [ (sequence, transcriptName)]
    val intmers = combined.map(tup => (IntMer.fromSequence(tup._1), tup._2)) // [ ( [Intmers] , transcriptName ) ]
    val kmers = intmers.flatMap(tup => { tup._1.map(imer => ((imer.longHash, imer.isOriginal, {if(imer.isOriginal) imer.toCanonicalString else imer.toAntiCanonicalString}), tup._2) )} ) //[ ((hash, orig), transcriptName) ]
    val partialMap = kmers.groupBy(_._1) // Map[ (hash, orig) -> [ ((hash, orig), transcripts)] ]
    val idx = partialMap.mapValues(v => v.groupBy(w => w._2).mapValues(_.size.toLong)) // Map[ (hash, orig) -> Map[ transcriptName -> count ] ]
    idx
  }

  /**
   * Creates two Indexex, one using the Index function and the other by manually cutting kmers.
   *
   * @param sequences A list of Strings which describe the base sequences of the contig fragments we are 
   *                  emulating
   * @return A tuple of Indexes, one from the Index function and the other calculated manually
   */
  // Compute actual Index and the expected results for Index based on the input sequences
  def createIndex(sequences: Array[String]): (Map[(Long, Boolean), Map[String, Long]], Map[(Long, Boolean, String), Map[String, Long]]) = {
    // Name all the sequences in order
    val names = {for (i <- 0 to sequences.size) yield "seq" + i.toString}.toArray

    // Create a one contig fragment per sequence
    val frags = sc.parallelize( sequences.zip(names).map(c => createContigFragment(c._1, c._2)) )

    // Create an arbitrary set of transcripts (only required as an ancillary argument to Index)
    val transcripts = sc.parallelize( Seq(Transcript("one", Seq("one"), "gene1", true, Iterable[Exon](), Iterable[CDS](), Iterable[UTR]()),
      Transcript("two", Seq("two"), "gene1", true, Iterable[Exon](), Iterable[CDS](), Iterable[UTR]())) )

    // Compute Index
    val (imap, tmap) = Index(frags, transcripts)

    // Compute expected results
    val expected = expectedResults(sequences, names)

    (imap, expected)
  }

  /**
   * Compares the results of the calculated Index with the expected results
   *
   * @param recieved Map recieved from the Index call
   * @param expected Manually calculated Map to compare recieved results against
   * @return Returns True if the recieved Map is identical to the expected one, false otherwise.
   */
  def compareResults(recieved: Map[(Long, Boolean), Map[String, Long]], expected: Map[(Long, Boolean, String), Map[String, Long]]) : Boolean = {
    // Make a set of the kmers in Expected that matches the format of the kmers in Recieved:
    val formattedExpected = expected.keySet.map(k => (k._1, k._2))

    // Check Sizes
    val recievedSize = recieved.size 
    val expectedSize = expected.size
    val sizesMatch = recievedSize == expectedSize

    // Edge case if either Recieved or Expected is empty
    val sizeMsg = "Recieved Index of Size: " + recievedSize.toString + ", but expected size of " + expectedSize.toString
    if (recievedSize == 0 || expectedSize == 0) {
      println("Index size was zero")
      println(sizeMsg)
      return false
    }

    // Check if all kmers in expected are present in recieved 
    val kmersPresent = recieved.keySet == formattedExpected

    // Check if there are kmers missing
    val missingKmers = expected.keySet.filter(k => {
      val key = (k._1, k._2)
      val inRecieved = recieved.keySet.contains(key)
      !inRecieved
    })
    val missingMsg = if (missingKmers.size > 0) "Kmers missing from Index:\n" + missingKmers.map(k => k.toString + "\n").reduce(_+_) else "No missing kmers"

    // Check if there were kmers added
    val addedKmers = recieved.keySet.filter(k => {
      val inExpected = formattedExpected.contains(k)
      !inExpected
    })
    val addedMsg = if (addedKmers.size > 0) "Kmers added to Index (that shouldn't be present):\n" + addedKmers.map(k => k.toString + "\n").reduce(_+_) else "No added kmers"

    val kmersCorrect = sizesMatch && kmersPresent
    val kmerMsg = sizeMsg + "\n" + missingMsg + "\n" + addedMsg + "\n"

    // Check if counts on the transcripts are correct:
    val transcriptsCheck = expected.keySet.map(e => {
      // For each expected kmer, check if it exists in recieved
      val key = (e._1, e._2)
      val keyInRecieved = recieved.contains(key)
      if (keyInRecieved) {
        // If the kmer does exist, then check if the expected transcripts match the recieved ones
        val recMap = recieved(key)
        val expMap = expected(e)
        val transcriptsCorrect = expMap.keySet.map(t => {
          val transcriptInRecieved = recMap.contains(t)
          if (transcriptInRecieved) {
            // If the transcript exists, then check if the count is accurate
            val recCount = recMap(t)
            val expCount = expMap(t)
            if (recCount == expCount) {
              (true, "")
            }
            else {
              val transcriptsMsg = "Count for transcript=" + t.toString + " for kmer=" + e.toString + " was " + recCount.toString + ". It should be " + expCount.toString
              (false, transcriptsMsg)
            }
          }
          else {
            val transcriptsMsg = t.toString + " not part of Map for " + e.toString + ". It should be."
            (false, transcriptsMsg)
          }
          })
        transcriptsCorrect.map(t => (t._1, if(t._2.size > 0) t._2 + "\n" else "")).reduce((t1, t2) => (t1._1 && t2._1, t1._2 + t2._2))
      }
      else {
        val transcriptsMsg = e.toString + " not in Index. It should be."
        (false, transcriptsMsg)
      }
      })
    val (transcriptsCountsCorrect, transcriptsMsg) = transcriptsCheck.map(t => (t._1, if(t._2.size > 0) t._2 + "\n" else "")).reduce((t1, t2) => (t1._1 && t2._1, t1._2 + t2._2))

    val correct = kmersCorrect && transcriptsCountsCorrect
    val msg = "############ Kmer Statistics ############\n" + kmerMsg + "\n############ Transcript Statistics ############\n" + transcriptsMsg
    val fullDump = "Expected Results:\n" + expected.map(q => q.toString + "\n").reduce(_ + _) + "Recieved Results:\n" + recieved.map(q => q.toString + "\n").reduce(_ + _)
    if (correct) {
      true
    }
    else {
      println(msg)
      println(fullDump)
      false
    }
  }

  /**
   * For a given set of input sequences, tests if Index produces the correct output
   *
   * @param sequences A list of Strings which describe the base sequences of the contig fragments we are 
   *                  emulating
   * @return Returns True if the sequences are correctly indexed, false otherwise.
   */
  def testOfIndex(sequences: Array[String]): Boolean = {
    val (recieved, expected) = createIndex(sequences)
    val correct = compareResults(recieved, expected)
    correct
  }

  sparkTest("simple sequences correctly indexed") { 

    val seq = "ACACTGTGGGTACACTACGAGA"
    val testPassed = testOfIndex(Array(seq))
    assert(testPassed)
  }

  sparkTest("repetetive sequences correctly indexed") {
    // Repeats of AAAAAAAAAAAAAAAA
    val seq1_1 = "AAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAA"
    val seq1_2 = "AAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTT"
    val testPassed1 = testOfIndex(Array(seq1_1, seq1_2))
    assert(testPassed1)
  }

  sparkTest("random sequences correctly indexed") {
    // Randomized string of length 32
    val seq2_1 = "GGATCAAATACGGGGCCTGTGTGCTGCGACTA"
    val seq2_2 = "ACTAGGGCCTGCATGCGAACATCCTGAACGCC"
    val testPassed2 = testOfIndex(Array(seq2_1, seq2_2))
    assert(testPassed2)
  }

  sparkTest("kmers with Canonicality=False make it through") {
    val seq1_1 = "TTTTTTTTTTTTTTTT"
    val testPassed1 = testOfIndex(Array(seq1_1))
    assert(testPassed1)

    val seq2_1 = "TTTTTTTTTTTTTTTT"
    val seq2_2 = "AAAAAAAAAAAAAAAA"
    val testPassed2 = testOfIndex(Array(seq2_1, seq2_2))
    assert(testPassed2)
  }
  
}
