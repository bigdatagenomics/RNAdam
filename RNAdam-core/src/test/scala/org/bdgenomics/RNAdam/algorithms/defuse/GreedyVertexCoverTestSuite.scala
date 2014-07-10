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

import org.bdgenomics.adam.util.SparkFunSuite

class GreedyVertexCoverTestSuite extends SparkFunSuite {

  val coverAlgorithm = new GreedyVertexCover()

  sparkTest("a universe with one subset results in that subset being chosen") {

    val universe = sc.parallelize(Seq(1, 2, 3))
    val subsets = sc.parallelize(Seq(0L -> Set(1, 2, 3)))

    val assignments: Map[Int, Long] =
      coverAlgorithm.calculateSetCover[Int](universe, subsets).collect().toMap

    assert(assignments.size === 3)
    assert(assignments(1) === 0L)
    assert(assignments(2) === 0L)
    assert(assignments(3) === 0L)
  }

  sparkTest("a universe with three subsets, one of which is redundant, results in only two being chosen") {
    val universe = sc.parallelize(Seq(1, 2, 3, 4, 5, 6))
    val subsets = sc.parallelize(Seq(0L -> Set(1, 2, 3, 4), 1L -> Set(5, 6), 2L -> Set(6)))

    val assignments: Map[Int, Long] =
      coverAlgorithm.calculateSetCover[Int](universe, subsets).collect().toMap

    assert(assignments.size === 6)
    assert(assignments(1) === 0L)
    assert(assignments(2) === 0L)
    assert(assignments(3) === 0L)
    assert(assignments(4) === 0L)
    assert(assignments(5) === 1L)
    assert(assignments(6) === 1L)
  }

  sparkTest("overlapping subsets can invalidate each other") {
    val universe = sc.parallelize(Seq(1, 2, 3, 4, 5, 6))
    val subsets = sc.parallelize(Seq(0L -> Set(1, 2, 3, 4), 1L -> Set(3, 4, 5), 2L -> Set(6)))

    val assignments: Map[Int, Long] =
      coverAlgorithm.calculateSetCover[Int](universe, subsets).collect().toMap

    assert(assignments.size === 6)
    assert(assignments(1) === 0L)
    assert(assignments(2) === 0L)
    assert(assignments(3) === 1L)
    assert(assignments(4) === 1L)
    assert(assignments(5) === 1L)
    assert(assignments(6) === 2L)
  }

}
