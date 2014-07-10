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
import org.apache.spark.graphx.{ EdgeTriplet, Graph, Pregel }
import org.apache.spark.rdd.RDD
import scala.reflect.ClassTag

object FindCliques extends Serializable {

  /**
   * Finds maximal cliques in a graph.
   *
   * @param graph Graph to find cliques on.
   * @return An RDD containing sequences of nodes in cliques.
   *
   * @tparam ED Edge data.
   * @tparam VD Vertex data.
   */
  def findMaximalCliques[ED, VD](graph: Graph[VD, ED])(implicit argED: ClassTag[ED], argVD: ClassTag[VD]): RDD[Seq[VD]] = {
    def labelNeighbors(graph: Graph[VD, ED]): Graph[CliqueVertex[VD], ED] = {
      // wrap vertices
      val cvGraph = graph.mapVertices((vid, v) => CliqueVertex(v, Seq(), Seq(), vid))

      // pregel to get neighbors
      Pregel(cvGraph, NeighborMessage(), 1)((vid: Long,
        v: CliqueVertex[VD],
        msg: NeighborMessage) => {
        CliqueVertex(v.data, msg.vid, Seq(Clique(Set(vid))), vid)
      }, et => {
        Iterator((et.dstId, NeighborMessage(et.srcId)),
          (et.srcId, NeighborMessage(et.dstId)))
      }, NeighborMessage.merge(_, _))
    }

    def labelVerticesWithCliques(graph: Graph[CliqueVertex[VD], ED]): Graph[CliqueVertex[VD], ED] = {
      // pregel to label cliques
      Pregel(graph, CliqueMessage(), 1)((vid: Long,
        v: CliqueVertex[VD],
        msg: CliqueMessage) => {
        CliqueVertex(v.data, v.neighbors, msg.cliques ++ v.cliques, vid)
      }, et => {
        Iterator((et.dstId, CliqueMessage(et.srcId, et.srcAttr.neighbors)),
          (et.srcId, CliqueMessage(et.dstId, et.dstAttr.neighbors)))
      }, CliqueMessage.merge(_, _))
    }

    // generate graphs
    val neighborGraph = labelNeighbors(graph)
    val labeledGraph = labelVerticesWithCliques(neighborGraph)

    // group by cliques
    labeledGraph.mapVertices((_, v) => v.reduce())
      .vertices
      .flatMap(kv => kv._2.mapCliques)
      .groupByKey()
      .map(kv => kv._2.toSeq.distinct)
  }

}
