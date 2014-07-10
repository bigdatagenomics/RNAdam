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
  override def calculateSetCover[U, S](universe: RDD[U], subsets: RDD[(S, Set[U])])(implicit argU: ClassTag[U], argS: ClassTag[S]): RDD[(U, S)] = {

    var assignment: RDD[(U, Option[S])] = universe.map(u => (u, None))

    while (assignment.filter(!_._2.isDefined).count() != 0) {

      // Find the as-yet-uncovered elements of the universe.
      val uncovered: RDD[(U, Option[S])] = assignment.filter(!_._2.isDefined)

      // For each subset (indexed by its subsetId), count how many of the
      // uncovered elements of the universe are covered by the subset
      // Call this the 'subset-coverage value' for each subset

      def mapWithSubset(pair: (S, Set[U])): Seq[(U, S)] = pair._2.map(u => (u, pair._1)).toSeq
      val elementToSubsetMap: RDD[(U, S)] = subsets.flatMap(mapWithSubset)

      val uncoveredMap: RDD[(U, (S, Option[S]))] = elementToSubsetMap.join(uncovered)

      val subsetToUncoveredMap: RDD[(S, U)] = uncoveredMap.map(p => (p._2._1, p._1))

      val subsetToGroupedUncovered: RDD[(S, Iterable[U])] = subsetToUncoveredMap.groupByKey()

      def countIterable(pair: (S, Iterable[U])): (S, Int) = (pair._1, pair._2.size)
      val subsetCoverage: RDD[(S, Int)] = subsetToGroupedUncovered.map(countIterable)

      // Find the MAX subset-coverage value
      val maxCoverage: Int = subsetCoverage.map(_._2).max()

      // Find the subset whose subset-coverage value is equal to the max --
      // if there is more than one such subset, pick one at random
      def hasMaxCoverage(pair: (S, Int)): Boolean = pair._2 >= maxCoverage
      val maxSubset: S = subsetCoverage.filter(hasMaxCoverage).first()._1

      // these are the elements in the universe which are now covered
      // by our chosen maximal subset.
      val covered: RDD[(U, S)] = subsets.filter(_._1 == maxSubset).flatMap(mapWithSubset)

      def combineAssignments(assignment: (U, (Option[S], Option[S]))): (U, Option[S]) =
        assignment._2._2 match {
          case None                => (assignment._1, assignment._2._1)
          case Some(newAssignment) => (assignment._1, Some(newAssignment))
        }

      // Update the assignments, to include all the new elements.
      assignment = assignment.leftOuterJoin(covered).map(combineAssignments)
    }

    assignment.map(p => (p._1, p._2.get))
  }

}
