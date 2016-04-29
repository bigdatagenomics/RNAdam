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
package org.bdgenomics.rice

import org.bdgenomics.utils.instrumentation.Metrics

/**
 * Contains [[Timers]] that are used to instrument rice.
 */
private[rice] object Timers extends Metrics {

  // CLI
  val LoadingContigs = timer("Loading Contig Fragments")
  val LoadingTranscripts = timer("Loading Transcripts")
  val Indexing = timer("Indexing k-mers")
  val Saving = timer("Saving Index to disk")

  // Indexing
  val ContigsToGraph = timer("Converting to DeBruijnGraph")
  val VertexMapping = timer("Creating a map of kmers to counts")

  // Quantification
  val ExtractTranscriptLengths = timer("Extraction Transcript Lengths")
  val CountKmers = timer("Counting k-mers")
  val TareKmers = timer("Calibrate k-mer counts vs. GC Content")
  val CountEquivalenceClasses = timer("Map k-mers to Equivalence Classes")
  val NormalizingCounts = timer("Normalizing Input Counts")
  val InitializingEM = timer("Initializing EM Algorithm")
  val InitializingCounts = timer("Initializing Counts")
  val InitializingMu = timer("Initializing µ's")
  val RunningEMIter = timer("Running an Iteration of EM")
  val EStage = timer("E stage of EM")
  val MStage = timer("M stage of EM")
  val CalibratingForLength = timer("Calibrating vs. Transcript Length")
  val JoiningAgainstTranscripts = timer("Joining vs. Transcripts")
}
