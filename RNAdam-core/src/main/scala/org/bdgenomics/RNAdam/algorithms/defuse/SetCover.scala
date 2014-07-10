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

import scala.reflect.ClassTag

trait SetCover extends Serializable {

  /**
   * Calculates an answer to the 'set cover' problem (c.f. http://en.wikipedia.org/wiki/Set_cover_problem)
   * for a given total set ('universe') of elements.
   *
   * @param universe The complete set of elements to be covered.
   * @param subsets An (indexed) set of subsets of the elements in the universe
   * @return An assignment of each element of the universe to one of the indexed subsets.
   */
  def calculateSetCover[U](universe: RDD[U], subsets: RDD[(Long, Set[U])])(implicit argU: ClassTag[U]): RDD[(U, Long)]
}
