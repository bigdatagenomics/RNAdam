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
package org.bdgenomics.RNAdam.algorithms.defuse

import org.apache.spark.rdd.RDD
import org.bdgenomics.RNAdam.models.ReadPair
import org.bdgenomics.adam.rich.RichADAMRecord._

object FragmentLengthDistribution extends Serializable {

  /**
   * Calculates a fragment length distribution, and excludes outliers given an
   * alpha parameter.
   *
   * @param reads An RDD of ADAM reads.
   * @param alpha The top/bottom % of reads to exclude.
   * @param sampleBy Optional down-sampling parameter. Default is not to sample.
   * @return (l_{min}, l_{max}): Return the min and max length.
   */
  def findPercentiles(reads: RDD[ReadPair],
                      alpha: Double,
                      sampleBy: Option[Double] = None): (Long, Long) = {
    // downsample our reads if a sampling level is provided
    val sampledReads = sampleBy.fold(reads)(reads.sample(false, _))

    // calculate the insert distribution
    val insertDistribution = sampledReads.map(pair => {
      // get end of first read and start of second read
      val endFirst = pair.first.end
      val startSecond = Option(pair.second.getStart)
      assert(endFirst.isDefined && startSecond.isDefined,
        "Reads don't seem to be aligned...")

      // calculate insert size
      startSecond.get - endFirst.get
    }).countByValue()

    // calculate the number of reads we have
    val numSamples = insertDistribution.values.reduce(_ + _)

    // sort the insert size/count distribution
    val sampleSeq = insertDistribution.toSeq
      .sortBy(kv => kv._1)

    // calculate min/max counts by dropping until we've seen enough
    var (startCount, endCount) = (0L, 0L)
    val min = sampleSeq.dropWhile(kv => {
      startCount += kv._2
      kv._1 < numSamples * alpha
    }).head._1
    val max = sampleSeq
      .reverse
      .dropWhile(kv => {
        startCount += kv._2
        kv._1 < numSamples * alpha
      }).head._1

    // return
    (min, max)
  }
}
