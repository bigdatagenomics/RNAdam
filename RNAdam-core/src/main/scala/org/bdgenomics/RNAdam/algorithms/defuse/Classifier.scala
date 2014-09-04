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
import org.bdgenomics.RNAdam.models.ReadPair
import org.bdgenomics.formats.avro.ADAMRecord

object Classifier {

  def classify(records: RDD[ADAMRecord]): (RDD[ReadPair], RDD[ReadPair], RDD[ReadPair]) = {

    val r1: RDD[(String, ADAMRecord)] = records.keyBy(x => x.getReadName.toString)

    val groupedByReadName: RDD[(String, Seq[ADAMRecord])] = r1.groupBy(p => p._1).map {
      case (key: String, iter: Iterable[(String, ADAMRecord)]) => (key, iter.map(x => x._2).toSeq)
    }

    val readPairs: RDD[ReadPair] = groupedByReadName.flatMap {
      case (key: String, records: Seq[ADAMRecord]) =>
        findReadPairs(records)
    }

    def concordant: RDD[ReadPair] = readPairs.filter(x => sameTranscript(x))
    def spanning: RDD[ReadPair] = readPairs.filter(x => spanningTranscript(x))
    def split: RDD[ReadPair] = readPairs.filter(x => splitTranscript(x))
    (concordant, spanning, split)
  }

  private def findReadPairs(records: Seq[ADAMRecord]): Seq[ReadPair] = {

    val firstRecords: Seq[ADAMRecord] = records.filter(_.getFirstOfPair)
    val secondRecords: Seq[ADAMRecord] = records.filter(_.getSecondOfPair)

    firstRecords.flatMap {
      case first: ADAMRecord =>
        secondRecords.map {
          case second: ADAMRecord =>
            ReadPair(first, second)
        }
    }
  }

  private def hasTranscriptName(record: ADAMRecord): Boolean = {
    record.getContig() != null
  }

  private def sameTranscript(pair: ReadPair): Boolean = {

    if (!hasTranscriptName(pair.first) || !hasTranscriptName(pair.second))
      return false
    pair.first.getContig.getContigName.equals(pair.second.getContig().getContigName)
  }

  private def spanningTranscript(pair: ReadPair): Boolean = {

    if (!hasTranscriptName(pair.first) || !hasTranscriptName(pair.second))
      return false
    !pair.first.getContig.getContigName.equals(pair.second.getContig().getContigName)
  }

  private def splitTranscript(pair: ReadPair): Boolean = {
    hasTranscriptName(pair.first) && !hasTranscriptName(pair.second)
  }
}
