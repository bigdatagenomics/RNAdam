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
package org.bdgenomics.rice.algorithms.alignment

import net.fnothaft.ananas.models.CanonicalKmer
import org.bdgenomics.rice.models.KmerIndex

object SimpleAlignmentModel extends AlignmentModel {

  def likelihood(count: Long, readLength: Long, kmerSize: Long = 16L): Double = {
    return count.toDouble / (readLength - kmerSize + 1)
  }

  /**
   * Given a read, returns (transcript, likelihood) pairs. ASSUMES KMER LENGTH IS 16
   */
  def processRead(iter: Iterator[CanonicalKmer],
                  kmerIndex: KmerIndex): Map[String, Double] = {

    val ar = iter.toArray // So we can use this more than once

    val kLength = ar(0).kmerLength // kmerLength
    val readLength = ar.map(c => 1).reduce(_ + _) + kLength - 1 // Length of read

    ar.flatMap(c => kmerIndex.getTranscripts(c)) // Iter[ (TranscriptId, Count) ]
      .map(c => (c._1, 1L)) // Iter[ (TranscriptId, 1L)]
      .groupBy(_._1) // Map[ TranscriptId -> Array(TranscriptId, 1) ]
      .map(t => (t._1, likelihood(t._2.size, readLength, kLength))) // Map[ TranscriptId -> likelihood ]
  }
}