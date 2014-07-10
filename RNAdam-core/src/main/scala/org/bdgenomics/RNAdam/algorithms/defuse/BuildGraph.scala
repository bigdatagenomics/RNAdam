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

import org.apache.spark.graphx.{ Edge, Graph }
import org.apache.spark.rdd.RDD
import org.bdgenomics.RNAdam.algorithms.cliques.FindCliques
import org.bdgenomics.RNAdam.models.{ ApproximateFusionEvent, ReadPair }

object BuildGraph extends Serializable {

  def apply(spanningRecords: RDD[ReadPair],
            lMax: Long): RDD[(ApproximateFusionEvent, Seq[ReadPair])] = {
    val initGraph = buildGraph(spanningRecords, lMax)
    val trimmedGraph = trimGraph(initGraph)
    val cliques = FindCliques.findMaximalCliques(trimmedGraph)
    cliques.keyBy(s => {
      s.map(_.generateEvent(lMax))
        .reduce(_.intersection(_))
    })
  }

  /**
   * Trims edges frm the graph if the edges represent fusion events that don't overlap.
   *
   * @param graph Graph to trim.
   * @return Trimmed graph.
   */
  protected def trimGraph(graph: Graph[ReadPair, (ApproximateFusionEvent, ApproximateFusionEvent)]): Graph[ReadPair, (ApproximateFusionEvent, ApproximateFusionEvent)] = {
    graph.subgraph(et => {
      val fusionEvents = et.attr
      // do the two approximate fusion events overlap?
      fusionEvents._1.overlaps(fusionEvents._2)
    })
  }

  /**
   * Builds an unpruned graph out of read pairs.
   *
   * @param spanningRecords Read pairs that span fusion points.
   * @param lMax Max window length.
   * @return Returns an unpruned graph.
   */
  protected def buildGraph(spanningRecords: RDD[ReadPair],
                           lMax: Long): Graph[ReadPair, (ApproximateFusionEvent, ApproximateFusionEvent)] = {

    // group records to split transcripts
    val groupedRecords = spanningRecords.groupBy(rp => {
      // get transcript names
      val firstTranscript = rp.first.getContig.getContigName.toString
      val secondTranscript = rp.second.getContig.getContigName.toString
      assert(firstTranscript != secondTranscript,
        "Spanning reads should be on different transcripts.")

      // return the two transcript names joined together
      firstTranscript + "," + secondTranscript
    })

    // generate edges
    val recordsAndEdges = groupedRecords.flatMap(kv => {
      val (_, readPairs) = kv

      // generate ids and the events
      val idsAndEvents = readPairs.map(r => (r.hash, r.generateEvent(lMax))).toArray

      // initialize
      var edges = List[Edge[(ApproximateFusionEvent, ApproximateFusionEvent)]]()

      // loop and pair
      (0 until idsAndEvents.size).foreach(i => {
        (i + 1 until idsAndEvents.size).foreach(j => {
          val src = idsAndEvents(i)
          val dst = idsAndEvents(j)

          edges = Edge(src._1, dst._1, (src._2, dst._2)) :: edges
        })
      })

      // return
      edges
    })

    // build graph
    Graph(spanningRecords.keyBy(_.hash), recordsAndEdges)
  }
}
