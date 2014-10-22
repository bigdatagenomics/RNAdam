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

package org.bdgenomics.RNAdam.algorithms.quantification

import org.bdgenomics.adam.util.SparkFunSuite
import org.apache.spark.Logging
import org.apache.spark.rdd.RDD
import org.apache.spark.SparkContext._
import org.bdgenomics.adam.models._
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rdd.read.AlignmentRecordContext._
import org.bdgenomics.formats.avro.AlignmentRecord
import org.bdgenomics.RNAdam.algorithms.quantification.Quantify._

class QuantificationTestSuite extends SparkFunSuite {

  def dummyTranscript(id: String): Transcript = {
    return new Transcript(id,
      Seq("test"),
      "Gene1",
      true,
      List[Exon](),
      List[CDS](),
      List[UTR]());
  }

  sparkTest("Dummy Transcript correctly initialized:") {
    var t = dummyTranscript("t1");
    assert(t.id == "t1")
    assert(t.strand == true)
  }

  sparkTest("transcripts correctly matched with coverage") {
    var s1: Double = 1
    var s2: Double = 2
    var s3: Double = 3
    val rdd1 = sc.parallelize(Array(dummyTranscript("t1"), dummyTranscript("t2"), dummyTranscript("t3")))
    val rdd2 = sc.parallelize(Array(("t1", s1, Iterable[Long]()), ("t2", s2, Iterable[Long]()), ("t3", s3, Iterable[Long]())))
    val rdd3 = Quantify.joinTranscripts(rdd1, rdd2)
    val compare = rdd3.collect()
    for (x <- compare) {
      if (x._1.id == "t1") {
        assert(x._2 == s1)
      }
      if (x._1.id == "t2") {
        assert(x._2 == s2)
      }
      if (x._1.id == "t3") {
        assert(x._2 == s3)
      }

    }

  }

}

