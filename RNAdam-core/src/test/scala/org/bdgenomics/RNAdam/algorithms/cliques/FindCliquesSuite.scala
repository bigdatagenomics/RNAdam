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
package org.bdgenomics.RNAdam.algorithms.cliques

import org.apache.spark.SparkContext._
import org.apache.spark.graphx.{ Edge, Graph }
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.util.SparkFunSuite

class FindCliquesSuite extends SparkFunSuite {

  sparkTest("disconnected nodes are their own clique") {
    val graph = Graph(sc.parallelize(Seq((0L, 0),
      (1L, 1),
      (2L, 2))),
      sc.parallelize(Seq[Edge[Int]]()))

    val cliques = FindCliques.findMaximalCliques(graph)

    assert(cliques.count() === 3)
    cliques.collect.foreach(i => assert(i.length === 1))
  }

  sparkTest("gather small but disconnected cliques") {
    val graph = Graph(sc.parallelize(Seq((0L, 0),
      (1L, 1),
      (2L, 2),
      (3L, 3))),
      sc.parallelize(Seq[Edge[Int]](Edge(0L, 1L, 0),
        Edge(2L, 3L, 0))))

    val cliques = FindCliques.findMaximalCliques(graph)

    assert(cliques.count() === 2)
    cliques.collect.foreach(i => assert(i.length === 2))
  }

  sparkTest("gather small, overlapping cliques") {
    val graph = Graph(sc.parallelize(Seq((0L, 0),
      (1L, 1),
      (2L, 2),
      (3L, 3))),
      sc.parallelize(Seq[Edge[Int]](Edge(0L, 1L, 0),
        Edge(0L, 2L, 0),
        Edge(1L, 2L, 0),
        Edge(2L, 3L, 0))))

    val cliques = FindCliques.findMaximalCliques(graph)

    assert(cliques.count() === 2)
    val sizes = cliques.map(_.size)
      .countByValue
    assert(sizes(3) === 1)
    assert(sizes(2) === 1)
  }

}
