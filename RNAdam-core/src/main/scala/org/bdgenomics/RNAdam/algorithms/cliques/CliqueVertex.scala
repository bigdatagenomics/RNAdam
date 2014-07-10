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

private[cliques] case class CliqueVertex[VD](data: VD,
                                             neighbors: Seq[Long],
                                             cliques: Seq[Clique],
                                             vid: Long) {

  override def toString(): String = {
    "value: " + data.toString + " neighbors: " + neighbors.toString + " cliques: " + cliques
  }

  /**
   * Reduces down the set of candidate cliques to the set of actual cliques.
   *
   * @return Returns a reduced clique vertex.
   */
  def reduce(): CliqueVertex[VD] = {
    val reducedCliques = cliques.map(c => c.nodes.filter(v => neighbors.contains(v) || v == vid))

    val newCliques = reducedCliques.filter(c => reducedCliques.forall(v => !(c.subsetOf(v) &&
      !c.equals(v))))
      .map(Clique(_))

    CliqueVertex(data, neighbors, newCliques, vid)
  }

  /**
   * Used for multi-mapping a vertex to all the cliques it belongs in.
   *
   * @return A seq of (clique, vertex) mappings.
   */
  def mapCliques(): Seq[(Clique, VD)] = cliques.map(c => (c, data))
}
