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
import org.apache.spark.SparkContext._

import scala.reflect.ClassTag

class GreedyVertexCover extends SetCover {

  /**
   * Calculates an answer to the 'set cover' problem (c.f. http://en.wikipedia.org/wiki/Set_cover_problem)
   * for a given total set ('universe') of elements.
   *
   * @param universe The complete set of elements to be covered.
   * @param subsets An (indexed) set of subsets of the elements in the universe
   * @return A set of indices -- taking the union of all subsets indexed by the result should
   *         recover the universe
   */
  override def calculateSetCover[U](universe: RDD[U], subsets: RDD[(Long, Set[U])])(implicit argU: ClassTag[U]): RDD[(U, Long)] = {

    var assignment: RDD[(U, Long)] = universe.map(u => (u, -1L))

    while (assignment.filter(_._2 == -1).count() != 0) {

      // Find the as-yet-uncovered elements of the universe.
      val uncovered: RDD[(U, Long)] = assignment.filter(_._2 == -1L)

      // For each subset (indexed by its subsetId), count how many of the
      // uncovered elements of the universe are covered by the subset
      // Call this the 'subset-coverage value' for each subset
      val elementToSubsetMap: RDD[(U, Long)] = subsets.flatMap {
        case (id: Long, subset: Set[U]) => subset.map(u => (u, id)).toSeq
      }

      val uncoveredMap: RDD[(U, (Long, Long))] = elementToSubsetMap.join(uncovered)

      val subsetToUncoveredMap: RDD[(Long, U)] = uncoveredMap.map(p => (p._2._1, p._1))

      val subsetToGroupedUncovered: RDD[(Long, Iterable[U])] = subsetToUncoveredMap.groupByKey()

      val subsetCoverage: RDD[(Long, Int)] = subsetToGroupedUncovered.map {
        case (subsetId: Long, uncoveredElements: Iterable[U]) => (subsetId, uncoveredElements.size)
      }

      // Find the MAX subset-coverage value
      val maxCoverage: Int = subsetCoverage.map(_._2).max()

      // Find the subset whose subset-coverage value is equal to the max --
      // if there is more than one such subset, pick one at random
      val maxSubsetId: Long = subsetCoverage.filter {
        case (subsetId: Long, coverageSize: Int) => coverageSize >= maxCoverage
      }.first()._1

      // these are the elements in the universe which are now covered
      // by our chosen maximal subset.
      val covered: RDD[(U, Long)] = subsets.filter(_._1 == maxSubsetId).flatMap {
        case (id: Long, subset: Set[U]) => subset.map(u => (u, maxSubsetId))
      }

      def combineAssignments(assignment: (U, (Long, Option[Long]))): (U, Long) =
        assignment._2._2 match {
          case None => (assignment._1, assignment._2._1)
          case Some(newAssignment) => (assignment._1, newAssignment)
        }

      // Update the assignments, to include all the new elements.
      assignment = assignment.leftOuterJoin(covered).map(combineAssignments)
    }

    assignment
  }

}
