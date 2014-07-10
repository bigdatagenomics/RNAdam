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

private[cliques] object NeighborMessage {

  /**
   * Creates an empty neighbor message. Used for initial pregel step.
   *
   * @return Uninitialized neighbor message.
   */
  def apply(): NeighborMessage = NeighborMessage(Seq())

  /**
   * Returns a basic neighbor message.
   *
   * @param vid Vertex ID.
   * @return Returns a node specific neighbor message.
   */
  def apply(vid: Long): NeighborMessage = NeighborMessage(Seq(vid))

  /**
   * Merges two neighbor messages.
   *
   * @param n1 Message to merge.
   * @param n2 Message to merge.
   * @return Merged message.
   */
  def merge(n1: NeighborMessage, n2: NeighborMessage): NeighborMessage = {
    new NeighborMessage(n1.vid ++ n2.vid)
  }
}

private[cliques] case class NeighborMessage(vid: Seq[Long]) {
}
