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

import org.apache.spark.rdd.RDD
import org.bdgenomics.rice.utils.RiceFunSuite
import org.bdgenomics.rice.utils.TranscriptGenerator
import scala.math.{ abs, exp, log, max, min }
import scala.util.Random

class TareSuite extends RiceFunSuite {

  def fpEquals(a: Double, b: Double, eps: Double = 1e-6): Boolean = {
    val passed = abs(a - b) <= eps
    if (!passed) {
      println("|" + a + " - " + b + "| = " + abs(a - b) + "> " + eps)
    }
    passed
  }

  test("can't process illegal k-mers") {
    intercept[AssertionError] {
      Tare.kmerToDinucFeature("AN", 1L)
    }
    intercept[AssertionError] {
      Tare.kmerToDinucFeature("A", 10L)
    }
    intercept[AssertionError] {
      Tare.kmerToDinucFeature("ANTNC", 100L)
    }
  }

  test("chop a 2-mer into a feature") {
    val featureAA = Tare.kmerToDinucFeature("AA", 1000L)
    assert(fpEquals(featureAA.label, log(1000L)))
    assert(fpEquals(featureAA.features(0), 1.0))
    (1 to 15).foreach(i => assert(fpEquals(featureAA.features(i), 0.0)))

    val featureTT = Tare.kmerToDinucFeature("TT", 100L)
    assert(fpEquals(featureTT.label, log(100L)))
    assert(fpEquals(featureTT.features(15), 1.0))
    (0 to 14).foreach(i => assert(fpEquals(featureTT.features(i), 0.0)))
  }

  test("chop an 5-mer with a bad base into a feature") {
    val feature = Tare.kmerToDinucFeature("AANTT", 4321L)
    assert(fpEquals(feature.label, log(4321L)))
    assert(fpEquals(feature.features(0), 0.5))
    assert(fpEquals(feature.features(15), 0.5))
    (1 to 14).foreach(i => assert(fpEquals(feature.features(i), 0.0)))
  }

  sparkTest("generate biased kmers and try correcting their counts") {
    val sampleString = TranscriptGenerator.generateString(500, new Random(121212L))

    // generate kmers corresponding to log space gc bias curve of 2.0 + 1.0 (gc - 0.5)
    val kmerSamples = sampleString.sliding(15)
      .map(s => {
        // compute bias and use curve for count
        val gc = s.count(c => c == 'C' || c == 'G').toDouble / 15.0
        val count = (100.0 * exp(2.0 + 1.0 * (gc - 0.5))).toLong

        (s, count)
      }).toSeq

    // build rdd and fit sequence context curve
    val rdd = sc.parallelize(kmerSamples)
    val originalMax = rdd.map(kv => kv._2).reduce(_ max _)
    val originalMin = rdd.map(kv => kv._2).reduce(_ min _)
    val calibratedRdd = Tare.calibrateKmers(rdd)
      .cache()

    // calibration should reduce max and increase min
    val newMax = calibratedRdd.map(kv => kv._2).reduce(_ max _)
    val newMin = calibratedRdd.map(kv => kv._2).reduce(_ min _)

    assert(originalMax > newMax)
    assert(originalMin < newMin)
  }

  sparkTest("test of Tare.calibrateTxLenBias for 4 hand-picked values") {
    val itLong: Iterable[Long] = Seq(1L, 2L, 3L)
    val muHat: RDD[(String, Double, Iterable[Long])] = sc.parallelize(Seq(
      ("a", 0.28, itLong),
      ("b", 0.17, itLong),
      ("c", 0.31, itLong),
      ("d", 0.24, itLong)))
    val tLen: scala.collection.Map[String, Long] = (new scala.collection.immutable.HashMap[String, Long]()).+(
      ("a", 28L),
      ("b", 17L),
      ("c", 31L),
      ("d", 24L))
    val calMuHat: RDD[(String, Double, Iterable[Long])] = Tare.calibrateTxLenBias(muHat, tLen)
    assert(calMuHat.count() === 4)
    val a: Double = calMuHat.filter((x => x._1 == "a")).first()._2
    val b: Double = calMuHat.filter((x => x._1 == "b")).first()._2
    val c: Double = calMuHat.filter((x => x._1 == "c")).first()._2
    val d: Double = calMuHat.filter((x => x._1 == "d")).first()._2
    assert(fpEquals(a, 0.25))
    assert(fpEquals(b, 0.25))
    assert(fpEquals(c, 0.25))
    assert(fpEquals(d, 0.25))
  }

  sparkTest("randomized test of Tare.calibrateTxLenBias for small data size") {
    test(10)
  }

  sparkTest("randomize test of Tare.calibrateTxLenBias for larger data size") {
    test(10000)
  }

  def test(dataSize: Int) = {
    // This iterable of longs is just a placeholder.
    val itLong: Iterable[Long] = Seq(1L, 2L, 3L)

    // Randomly generate a sequence of numbers and use them to assemble data.
    // The seed value was chosen at 11:34 AM on Feb 6
    val rand: Random = new Random(113402062015L)
    val r: Seq[(String, Long)] = (0 to dataSize).map(i => (i.toString, 1L + rand.nextInt(10)))
    val sum: Double = r.map(x => x._2).reduce((x, y) => x + y).toDouble
    val tLen: scala.collection.Map[String, Long] = (new scala.collection.immutable.HashMap()).++(r)
    val muHat: RDD[(String, Double, Iterable[Long])] = sc.parallelize(r).map(x => (x._1, x._2 / sum, itLong))

    // Runs Tare.calibrateTxLenBias()
    val calMuHat: RDD[(String, Double, Iterable[Long])] = Tare.calibrateTxLenBias(muHat, tLen)

    // Checks for correct answer. Since all variation in abundance
    // is due to transcript length, all transcript should have equal
    // calibrated abundance.
    calMuHat.collect().foreach(x => assert(fpEquals(x._2, 1.0 / (dataSize + 1))))
  }

}
