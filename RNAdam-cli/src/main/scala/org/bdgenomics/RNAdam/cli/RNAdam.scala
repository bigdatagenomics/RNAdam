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
package org.bdgenomics.RNAdam.cli

import org.apache.hadoop.mapreduce.Job
import org.apache.spark.rdd.RDD
import org.apache.spark.{ SparkContext, Logging }
import org.kohsuke.args4j.{ Option => option, Argument }
import org.bdgenomics.formats.avro.{
  ADAMRecord,
  ADAMNucleotideContigFragment
}
import org.bdgenomics.adam.cli.{
  ADAMSparkCommand,
  ADAMCommandCompanion,
  ParquetArgs,
  SparkArgs,
  Args4j,
  Args4jBase
}
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.RNAdam.algorithms.defuse.Defuse
import org.bdgenomics.RNAdam.models.FusionEvent

object RNAdam extends ADAMCommandCompanion {

  val commandName = "RNAdam"
  val commandDescription = "Predict fusion transcripts using ADAM"

  def apply(args: Array[String]) = {
    new RNAdam(Args4j[RNAdamArgs](args))
  }
}

class RNAdamArgs extends Args4jBase with ParquetArgs with SparkArgs {
  @Argument(metaVar = "READS", required = true, usage = "ADAM read-oriented data", index = 0)
  var readInput: String = _

  @option(name = "-debug", usage = "If set, prints a higher level of debug output.")
  var debug = false

  @option(name = "-α", usage = "Alpha parameter for reads to discard for insert size calculation")
  var alpha = 0.0
}

class RNAdam(protected val args: RNAdamArgs) extends ADAMSparkCommand[RNAdamArgs] with Logging {

  // companion object to this class - needed for ADAMCommand framework
  val companion = RNAdam

  /**
   * Main method. Implements body of variant caller. SparkContext and Hadoop Job are provided
   * by the ADAMSparkCommand shell.
   *
   * @param sc SparkContext for RDDs.
   * @param job Hadoop Job container for file I/O.
   */
  def run(sc: SparkContext, job: Job) {

    log.info("Starting assembler...")

    log.info("Loading reads in from " + args.readInput)
    // load in reads from ADAM file
    val reads: RDD[ADAMRecord] = sc.adamLoad(args.readInput)

    // run defuse
    val fusions = Defuse.run(reads, args.alpha)

    // print fusions
    fusions.foreach(f => println("Fusion event at:" + f.start + "," + f.end))
  }
}
