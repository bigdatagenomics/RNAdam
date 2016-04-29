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
import org.bdgenomics.formats.avro.AlignmentRecord
import org.bdgenomics.rice.algorithms.alignment.AlignmentModel
import org.bdgenomics.rice.models.KmerIndex
import net.fnothaft.ananas.models.IntMer

object Mapper extends Serializable with Logging {

  def apply(reads: RDD[AlignmentRecord],
            kmerIndex: KmerIndex,
            likelihoodModel: AlignmentModel): RDD[(Long, Map[String, Double])] = {

    reads.zipWithUniqueId() // RDD[ read, readID ]
      .map(r => (r._2, likelihoodModel.processRead(IntMer.fromSequence(r._1.sequence).toIterator, kmerIndex))) // RDD[ readID, Map[color -> likelihood] ] 
  }
}
