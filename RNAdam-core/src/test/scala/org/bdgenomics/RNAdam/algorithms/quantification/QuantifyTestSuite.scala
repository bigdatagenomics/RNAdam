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

import scala.math
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.util.SparkFunSuite

class QuantifyTestSuite extends SparkFunSuite {

  sparkTest("test of mapKmersToClasses") {
    val kmerToEquivalenceClass: RDD[(String, Long)] = sc.parallelize(Seq(("a", 2),
      ("b", 3),
      ("c", 2),
      ("d", 1),
      ("e", 3)))
    val kmerCounts: RDD[(String, Long)] = sc.parallelize(Seq(("d", 80), ("a", 25), ("c", 35), ("b", 37), ("e", 38)))
    val classCounts: RDD[(Long, Long)] = Quantify.mapKmersToClasses(kmerCounts, kmerToEquivalenceClass)
    assert(classCounts.count() === 3)
    assert(classCounts.filter((x: (Long, Long)) => x._1 == 1).first() === (1, 80))
    assert(classCounts.filter((x: (Long, Long)) => x._1 == 2).first() === (2, 60))
    assert(classCounts.filter((x: (Long, Long)) => x._1 == 3).first() === (3, 75))
  }

  sparkTest("test of initializeEM") {
    val equivalenceClassCounts: RDD[(Long, Long)] = sc.parallelize(
      Seq((1, 45),
        (2, 52),
        (3, 49)))
    val equivalenceClassToTranscript: RDD[(Long, Iterable[String])] = sc.parallelize(
      Seq((2, Seq("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m")),
        (3, Seq("a", "b", "c", "d", "e", "f", "g")),
        (1, Seq("a", "b", "c", "d", "e"))))
    val result: RDD[(Long, Iterable[(String, Double)])] = Quantify.initializeEM(equivalenceClassCounts, equivalenceClassToTranscript)
    assert(result.count() === 3)
    val ec1p: RDD[(Long, Iterable[(String, Double)])] = result.filter((x: (Long, Iterable[(String, Double)])) => x._1 == 1)
    assert(ec1p.count() === 1)
    val ec1: (Long, Iterable[(String, Double)]) = ec1p.first()
    assert(ec1._1 === 1)
    assert(ec1._2.size === 5)
    assert(ec1._2.forall((x: (String, Double)) => {
      math.abs(x._2 - 9.0) < 1e-3
    }))
    val ec2p: RDD[(Long, Iterable[(String, Double)])] = result.filter((x: (Long, Iterable[(String, Double)])) => x._1 == 2)
    assert(ec2p.count() === 1)
    val ec2: (Long, Iterable[(String, Double)]) = ec2p.first()
    assert(ec2._1 === 2)
    assert(ec2._2.size === 13)
    assert(ec2._2.forall((x: (String, Double)) => {
      math.abs(x._2 - 4.0) < 1e-3
    }))
    val ec3p: RDD[(Long, Iterable[(String, Double)])] = result.filter((x: (Long, Iterable[(String, Double)])) => x._1 == 3)
    assert(ec3p.count() === 1)
    val ec3: (Long, Iterable[(String, Double)]) = ec3p.first()
    assert(ec3._1 === 3)
    assert(ec3._2.size === 7)
    assert(ec3._2.forall((x: (String, Double)) => {
      math.abs(x._2 - 7.0) < 1e-3
    }))
  }

}
