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
package org.bdgenomics.RNAdam.models

import org.bdgenomics.adam.models.ReferenceRegion
import org.bdgenomics.adam.rich.RichADAMRecord._
import org.bdgenomics.adam.rich.RichADAMRecord
import org.bdgenomics.formats.avro.ADAMRecord

case class ReadPair(first: ADAMRecord,
                    second: ADAMRecord) {

  /**
   * Generates a hash on this read pair.
   *
   * @return Returns a hash on this read pair's name and mapping site.
   */
  def hash(): Long = {
    def lString(s: String): Long = {
      // large prime
      var h = 1125899906842597L

      // loop and update
      s.foreach(c => h = 31 * h + c.toLong)

      // return
      h
    }

    // hash on read name, and transcript names
    lString(first.getReadName.toString + first.getContig.getContigName.toString +
      second.getContig.getContigName.toString)
  }

  /**
   * Generates an approximate fusion event from a read pair.
   *
   * @param l The length of the event window.
   * @return An approximate fusion event.
   *
   * @throws AssertionError Throws an assertion error if either read in the pair
   * is not aligned.
   */
  def generateEvent(l: Long): ApproximateFusionEvent = {
    val firstEnd = first.end
    val secondStart = Option(second.getStart)
    assert(firstEnd.isDefined && secondStart.isDefined,
      "Both reads in pair must be aligned.")

    // create fusion event
    ApproximateFusionEvent(ReferenceRegion(first.getContig.getContigName.toString,
      firstEnd.get, firstEnd.get + l),
      ReferenceRegion(second.getContig.getContigName.toString,
        secondStart.get - l, secondStart.get))
  }
}
