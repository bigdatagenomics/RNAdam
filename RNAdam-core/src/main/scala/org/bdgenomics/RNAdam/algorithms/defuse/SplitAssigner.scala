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
import org.bdgenomics.RNAdam.models.{ ApproximateFusionEvent, ReadPair }
import org.bdgenomics.adam.models.{ SequenceDictionary, ReferenceMapping, ReferenceRegion }
import org.bdgenomics.adam.rdd.RegionJoin
import org.bdgenomics.adam.rich.ReferenceMappingContext.ADAMRecordReferenceMapping
import org.bdgenomics.formats.avro.ADAMRecord
import scala.reflect._

object SplitAssigner {
  case object ReferenceRegionApproximateFusionEventReferenceMapping extends ReferenceMapping[(ReferenceRegion, ApproximateFusionEvent)] {
    override def getReferenceName(value: (ReferenceRegion, ApproximateFusionEvent)): String = value._1.referenceName

    override def getReferenceRegion(value: (ReferenceRegion, ApproximateFusionEvent)): ReferenceRegion = value._1
  }

  case object ADAMRecordReadPairReferenceMapping extends ReferenceMapping[(ADAMRecord, ReadPair)] {
    override def getReferenceName(value: (ADAMRecord, ReadPair)): String = ADAMRecordReferenceMapping.getReferenceName(value._1)

    override def getReferenceRegion(value: (ADAMRecord, ReadPair)): ReferenceRegion = ADAMRecordReferenceMapping.getReferenceRegion(value._1)
  }

  private[this] def referenceRegionsForEvents(events: RDD[ApproximateFusionEvent], lmin: Long, lmax: Long): RDD[(ReferenceRegion, ApproximateFusionEvent)] = {
    events.flatMap(afe => Seq((ReferenceRegion(afe.start.referenceName, afe.start.start - lmax, afe.start.end + lmax), afe),
      (ReferenceRegion(afe.end.referenceName, afe.end.start - lmax, afe.end.end + lmax), afe)))
  }

  def assignSplitsToFusions(events: RDD[ApproximateFusionEvent],
                            records: RDD[ReadPair],
                            seqDict: SequenceDictionary,
                            lmin: Long,
                            lmax: Long): RDD[(ApproximateFusionEvent, ReadPair)] = {
    val referenceRegions = referenceRegionsForEvents(events, lmin, lmax)
    val flattenedRecords = records.flatMap(rp => Seq(
      if (rp.first.getReadMapped) Some((rp.first, rp)) else None,
      if (rp.second.getReadMapped) Some((rp.second, rp)) else None)
      .flatten)
    RegionJoin.partitionAndJoin(events.sparkContext, seqDict, referenceRegions, flattenedRecords)(
      ReferenceRegionApproximateFusionEventReferenceMapping,
      ADAMRecordReadPairReferenceMapping,
      classTag[(ReferenceRegion, ApproximateFusionEvent)],
      classTag[(ADAMRecord, ReadPair)])
      .map {
        case ((_, afe), (_, rp)) => (afe, rp)
      }
  }
}
