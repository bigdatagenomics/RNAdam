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
package org.bdgenomics.rice.algorithms

import org.apache.spark.Logging
import org.apache.spark.rdd.RDD
import org.apache.spark.SparkContext._
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.models.Transcript
import org.bdgenomics.formats.avro.NucleotideContigFragment
import org.bdgenomics.rice.Timers._
import org.apache.spark.graphx.Graph
import net.fnothaft.ananas.debruijn.{ ColoredDeBruijnGraph, ColoredKmerVertex }
import net.fnothaft.ananas.models.ContigFragment

object Index extends Serializable with Logging {

  /**
   * Computes an index, given a set of Nucleotide Contig Fragments. An index provides a mapping 
   * from kmers to the transcripts they are found in 
   *
   * @param contigFragments An RDD containing contigFragments.
   * @param transcripts An RDD containing transcripts
   * @return Returns two mappings representing the Index of a sequence.
   */
  def apply(contigFragments: RDD[ContigFragment], transcripts: RDD[Transcript]): (Map[(Long, Boolean), Map[String, Long]], Map[String, Transcript]) = {

    val graph = createGraph(contigFragments)
    val vertexMapping = computeVertexMapping(graph)
    val transcriptMapping = computeTranscriptMapping(transcripts)
    (vertexMapping, transcriptMapping)
  }

  /**
   * Creates a colored de bruijn graph, given a set of Nucleotide Contig Fragments.
   *
   * @param contigFragments An RDD containing contigFragments.
   * @return Returns a Graph representing a colored De Bruijn graph of kmers
   */
  def createGraph(contigFragments: RDD[ContigFragment]): Graph[ColoredKmerVertex, Unit] = {

    ContigsToGraph.time {
      ColoredDeBruijnGraph.buildFromFragments(contigFragments)
    }
  }

  /**
   * Creates a mapping between kmers and the set of transcripts they appear in.
   *
   * @param graph A colored de bruijn graph representing kmers read from transcripts
   * @return Returns a Mapping from kmers to transcripts and the abundance of kmers in the transcripts
   */
  def computeVertexMapping(graph: Graph[ColoredKmerVertex, Unit]): Map[(Long, Boolean), Map[String, Long]] = {

    VertexMapping.time {

      // Convert each kmer in the graph into the tuple ( (kmerID, canonicality, transcriptItCameFrom) , 1 )
      val kmers: RDD[((Long, Boolean, String), Long)] = graph.vertices                        
        .flatMap(v => {
          // Grab forward kmers if they exist
          val forward = if (v._2.forwardTerminals.nonEmpty || v._2.forwardStronglyConnected.nonEmpty) {
            v._2.forwardTerminals.toList.map(t => ((v._1, true, t._1), 1L)) ++
              v._2.forwardStronglyConnected.toList.map(t => ((v._1, true, t._1._1), 1L))
          } else {
            Seq.empty
          }
          // Grab reverse kmers if they exist
          val reverse = if (v._2.reverseTerminals.nonEmpty || v._2.reverseStronglyConnected.nonEmpty) {
            v._2.reverseTerminals.toList.map(t => ((v._1, false, t._1), 1L)) ++
              v._2.reverseStronglyConnected.toList.map(t => ((v._1, false, t._1._1), 1L))
          } else {
            Seq.empty
          }

          // Combine foward and reverse
          forward ++ reverse
        }) 

      // Count the number of instances of each (kmerID, canonicality, transcriptItCameFrom) via reducing.
      val counts: RDD[((Long, Boolean, String), Long)] = kmers.reduceByKey(_ + _)

      // Reorder the tuples into ( (kmerID, canonicality) , (transcriptItCameFrom, count) )
      val reordered: RDD[((Long, Boolean), (String, Long))] = counts.map( q => {
          val ((hash, strand, transcript), number) = q
          ((hash, strand), (transcript, number))
        })

      // Collect all (kmerID, canonicality) by putting all associated (transcriptItCameFrom, count) into a map
      val collected: RDD[((Long, Boolean), Map[String, Long])] = reordered.combineByKey( v => 
          Map(v),
          (c, v) => c + v,
          (c1, c2) => c1 ++ c2 
        ) 

      // Convert the resulting rdd into a map from (kmerID, canonicality) => (transcriptItCameFrom, count)
      collected.collectAsMap().toMap 

    }
  }

  /**
   * Creates a Mapping between transcript IDs and Transcripts
   *
   * @param transcripts RDD of Transcripts
   * @return Returns a mapping between transcript IDs and transcript objects
   */
  def computeTranscriptMapping(transcripts: RDD[Transcript]): Map[String, Transcript] = {
    transcripts.map(t => (t.id, t)).collect().toMap
  }

}
