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
import org.bdgenomics.adam.models.{
  RecordGroupDictionary,
  ReferenceRegion,
  SequenceDictionary
}
import org.bdgenomics.adam.rdd.read.AlignmentRecordRDD
import org.bdgenomics.rice.models.{ Exon, Transcript }
import org.bdgenomics.rice.utils.riceFunSuite
import org.bdgenomics.rice.utils.{ ReadGenerator, TranscriptGenerator }
import scala.collection.Map
import scala.collection.immutable.HashMap
import scala.math.abs
import org.bdgenomics.adam.util.ReferenceFile
import org.bdgenomics.adam.util.{ TwoBitFile }
import org.bdgenomics.utils.io.{ ByteAccess, ByteArrayByteAccess }

class TestingTwoBitFile(byteAccess: ByteAccess) extends ReferenceFile with Serializable {
  // Test sequence, len = 24
  val testSeq = "CAATCCTTCGCCGCAGTGCA"
  override def extract(region: ReferenceRegion): String = {
    testSeq.substring(region.start.toInt, region.end.toInt)
  }
  override def sequences: SequenceDictionary = {
    SequenceDictionary.empty
  }
}

class QuantifySuite extends riceFunSuite {

  def fpEquals(a: Double, b: Double, eps: Double = 1e-6): Boolean = {
    val passed = abs(a - b) <= eps
    if (!passed) {
      println("|" + a + " - " + b + "| = " + abs(a - b) + "> " + eps)
    }
    passed
  }

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
      fpEquals(x._2, 9.0)
    }))
    val ec2p: RDD[(Long, Iterable[(String, Double)])] = result.filter((x: (Long, Iterable[(String, Double)])) => x._1 == 2)
    assert(ec2p.count() === 1)
    val ec2: (Long, Iterable[(String, Double)]) = ec2p.first()
    assert(ec2._1 === 2)
    assert(ec2._2.size === 13)
    assert(ec2._2.forall((x: (String, Double)) => {
      fpEquals(x._2, 4.0)
    }))
    val ec3p: RDD[(Long, Iterable[(String, Double)])] = result.filter((x: (Long, Iterable[(String, Double)])) => x._1 == 3)
    assert(ec3p.count() === 1)
    val ec3: (Long, Iterable[(String, Double)]) = ec3p.first()
    assert(ec3._1 === 3)
    assert(ec3._2.size === 7)
    assert(ec3._2.forall((x: (String, Double)) => {
      fpEquals(x._2, 7.0)
    }))
  }

  sparkTest("test of e") {
    // a, b, c, d are transcript names
    // 1, 2, 3, 4, 5, 6, 7 are equivalence class IDs
    val transcriptWeights: RDD[(String, Double, Iterable[Long])] = sc.parallelize(Seq(
      ("a", 2.0, Seq(1, 3, 5, 6)),
      ("b", 3.0, Seq(2, 4, 5)),
      ("c", 4.0, Seq(1, 2, 5, 6, 7)),
      ("d", 5.0, Seq(1, 2, 3))))
    val equivalenceClassAssignments: RDD[(Long, Iterable[(String, Double)])] = Quantify.e(transcriptWeights)
    assert(equivalenceClassAssignments.count() === 7)

    // variables to keep track of whether each transcript name appears exactly once, if it should appear at all
    var seen_a: Boolean = true
    var seen_b: Boolean = true
    var seen_c: Boolean = true
    var seen_d: Boolean = true

    // tests each record of the equivalenceClassAssignments RDD for correctness
    val record1RDD: RDD[(Long, Iterable[(String, Double)])] = equivalenceClassAssignments.filter((x: (Long, Iterable[(String, Double)])) => x._1 == 1)
    assert(record1RDD.count() === 1)
    val record1: (Long, Iterable[(String, Double)]) = record1RDD.first()
    assert(record1._1 === 1)
    assert(record1._2.toSeq.length === 3)
    seen_a = false
    seen_c = false
    seen_d = false
    record1._2.foreach((x: (String, Double)) => {
      if (x._1 == "a") {
        assert((!seen_a) && equalDouble(x._2, 2.0 / 11))
        seen_a = true
      } else if (x._1 == "c") {
        assert((!seen_c) && equalDouble(x._2, 4.0 / 11))
        seen_c = true
      } else if (x._1 == "d") {
        assert((!seen_d) && equalDouble(x._2, 5.0 / 11))
        seen_d = true
      } else {
        assert(false)
      }
    })

    val record2RDD: RDD[(Long, Iterable[(String, Double)])] = equivalenceClassAssignments.filter((x: (Long, Iterable[(String, Double)])) => x._1 == 2)
    assert(record2RDD.count() === 1)
    val record2: (Long, Iterable[(String, Double)]) = record2RDD.first()
    assert(record2._1 === 2)
    assert(record2._2.toSeq.length === 3)
    seen_b = false
    seen_c = false
    seen_d = false
    record2._2.foreach((x: (String, Double)) => {
      if (x._1 == "b") {
        assert((!seen_b) && equalDouble(x._2, 0.25))
        seen_b = true
      } else if (x._1 == "c") {
        assert((!seen_c) && equalDouble(x._2, 1.0 / 3))
        seen_c = true
      } else if (x._1 == "d") {
        assert((!seen_d) && equalDouble(x._2, 5.0 / 12))
        seen_d = true
      } else {
        assert(false)
      }
    })

    val record3RDD: RDD[(Long, Iterable[(String, Double)])] = equivalenceClassAssignments.filter((x: (Long, Iterable[(String, Double)])) => x._1 == 3)
    assert(record3RDD.count() === 1)
    val record3: (Long, Iterable[(String, Double)]) = record3RDD.first()
    assert(record3._1 === 3)
    assert(record3._2.toSeq.length === 2)
    seen_a = false
    seen_d = false
    record3._2.foreach((x: (String, Double)) => {
      if (x._1 == "a") {
        assert((!seen_a) && equalDouble(x._2, 2.0 / 7))
        seen_a = true
      } else if (x._1 == "d") {
        assert((!seen_d) && equalDouble(x._2, 5.0 / 7))
        seen_d = true
      } else {
        assert(false)
      }
    })

    val record4RDD: RDD[(Long, Iterable[(String, Double)])] = equivalenceClassAssignments.filter((x: (Long, Iterable[(String, Double)])) => x._1 == 4)
    assert(record4RDD.count() === 1)
    val record4: (Long, Iterable[(String, Double)]) = record4RDD.first()
    assert(record4._1 === 4)
    assert(record4._2.toSeq.length === 1)
    assert(record4._2.head._1 === "b")
    assert(equalDouble(record4._2.head._2, 1.0))

    val record5RDD: RDD[(Long, Iterable[(String, Double)])] = equivalenceClassAssignments.filter((x: (Long, Iterable[(String, Double)])) => x._1 == 5)
    assert(record5RDD.count() === 1)
    val record5: (Long, Iterable[(String, Double)]) = record5RDD.first()
    assert(record5._1 === 5)
    assert(record5._2.toSeq.length === 3)
    seen_a = false
    seen_b = false
    seen_c = false
    record5._2.foreach((x: (String, Double)) => {
      if (x._1 == "a") {
        assert((!seen_a) && equalDouble(x._2, 2.0 / 9))
        seen_a = true
      } else if (x._1 == "b") {
        assert((!seen_b) && equalDouble(x._2, 1.0 / 3))
        seen_b = true
      } else if (x._1 == "c") {
        assert((!seen_c) && equalDouble(x._2, 4.0 / 9))
        seen_c = true
      } else {
        assert(false)
      }
    })

    val record6RDD: RDD[(Long, Iterable[(String, Double)])] = equivalenceClassAssignments.filter((x: (Long, Iterable[(String, Double)])) => x._1 == 6)
    assert(record6RDD.count() === 1)
    val record6: (Long, Iterable[(String, Double)]) = record6RDD.first()
    assert(record6._1 === 6)
    assert(record6._2.toSeq.length === 2)
    seen_a = false
    seen_c = false
    record6._2.foreach((x: (String, Double)) => {
      if (x._1 == "a") {
        assert((!seen_a) && equalDouble(x._2, 1.0 / 3))
        seen_a = true
      } else if (x._1 == "c") {
        assert((!seen_c) && equalDouble(x._2, 2.0 / 3))
        seen_c = true
      } else {
        assert(false)
      }
    })

    val record7RDD: RDD[(Long, Iterable[(String, Double)])] = equivalenceClassAssignments.filter((x: (Long, Iterable[(String, Double)])) => x._1 == 7)
    assert(record7RDD.count() === 1)
    val record7: (Long, Iterable[(String, Double)]) = record7RDD.first()
    assert(record7._1 === 7)
    assert(record7._2.toSeq.length === 1)
    assert(record7._2.head._1 === "c")
    assert(equalDouble(record7._2.head._2, 1.0))

  }

  sparkTest("test of m") {
    val equivalenceClassAssignments: RDD[(Long, Iterable[(String, Double)])] = sc.parallelize(Seq(
      (1, Seq(("a", 0.6), ("c", 0.4))),
      (2, Seq(("b", 0.1), ("d", 0.5), ("a", 0.4))),
      (3, Seq(("a", 1.0))),
      (4, Seq(("c", 0.7), ("a", 0.3)))))
    val tLen: Map[String, Long] = (new HashMap()).+(
      ("a", 5),
      ("b", 6),
      ("c", 7),
      ("d", 3))
    val relNumKmersInEC: Map[Long, Double] = (new HashMap()).+(
      (1, 0.25),
      (2, 0.25),
      (3, 0.25),
      (4, 0.25))
    val transcriptWeights: RDD[(String, Double, Iterable[Long])] = Quantify.m(equivalenceClassAssignments, tLen, 3, relNumKmersInEC)
    assert(transcriptWeights.count() === 4)
    val a_record_RDD: RDD[(String, Double, Iterable[Long])] = transcriptWeights.filter((x: (String, Double, Iterable[Long])) => {
      x._1 == "a"
    })
    assert(a_record_RDD.count() === 1)
    val a_record: (String, Double, Iterable[Long]) = a_record_RDD.first()
    assert(a_record._1 === "a")
    assert(equalDouble(a_record._2, 460.0 / 907))
    val a_record_eq_cls: Seq[Long] = a_record._3.toSeq
    assert(a_record_eq_cls.length === 4)
    assert(a_record_eq_cls.filter((x: Long) => x == 1).length === 1)
    assert(a_record_eq_cls.filter((x: Long) => x == 2).length === 1)
    assert(a_record_eq_cls.filter((x: Long) => x == 3).length === 1)
    assert(a_record_eq_cls.filter((x: Long) => x == 4).length === 1)

    val b_record_RDD: RDD[(String, Double, Iterable[Long])] = transcriptWeights.filter((x: (String, Double, Iterable[Long])) => {
      x._1 == "b"
    })
    assert(b_record_RDD.count() === 1)
    val b_record: (String, Double, Iterable[Long]) = b_record_RDD.first()
    assert(b_record._1 === "b")
    assert(equalDouble(b_record._2, 15.0 / 907))
    val b_record_eq_cls: Seq[Long] = b_record._3.toSeq
    assert(b_record_eq_cls.length === 1)
    assert(b_record_eq_cls.filter((x: Long) => x == 1).length === 0)
    assert(b_record_eq_cls.filter((x: Long) => x == 2).length === 1)
    assert(b_record_eq_cls.filter((x: Long) => x == 3).length === 0)
    assert(b_record_eq_cls.filter((x: Long) => x == 4).length === 0)

    val c_record_RDD: RDD[(String, Double, Iterable[Long])] = transcriptWeights.filter((x: (String, Double, Iterable[Long])) => {
      x._1 == "c"
    })
    assert(c_record_RDD.count() === 1)
    val c_record: (String, Double, Iterable[Long]) = c_record_RDD.first()
    assert(c_record._1 === "c")
    assert(equalDouble(c_record._2, 132.0 / 907))
    val c_record_eq_cls: Seq[Long] = c_record._3.toSeq
    assert(c_record_eq_cls.length === 2)
    assert(c_record_eq_cls.filter((x: Long) => x == 1).length === 1)
    assert(c_record_eq_cls.filter((x: Long) => x == 2).length === 0)
    assert(c_record_eq_cls.filter((x: Long) => x == 3).length === 0)
    assert(c_record_eq_cls.filter((x: Long) => x == 4).length === 1)

    val d_record_RDD: RDD[(String, Double, Iterable[Long])] = transcriptWeights.filter((x: (String, Double, Iterable[Long])) => {
      x._1 == "d"
    })
    assert(d_record_RDD.count() === 1)
    val d_record: (String, Double, Iterable[Long]) = d_record_RDD.first()
    assert(d_record._1 === "d")
    assert(equalDouble(d_record._2, 300.0 / 907))
    val d_record_eq_cls: Seq[Long] = d_record._3.toSeq
    assert(d_record_eq_cls.length === 1)
    assert(d_record_eq_cls.filter((x: Long) => x == 1).length === 0)
    assert(d_record_eq_cls.filter((x: Long) => x == 2).length === 1)
    assert(d_record_eq_cls.filter((x: Long) => x == 3).length === 0)
    assert(d_record_eq_cls.filter((x: Long) => x == 4).length === 0)
  }

  def equalDouble(a: Double, b: Double): Boolean = {
    math.abs(a - b) < 1e-3
  }

  sparkTest("extract lengths from transcripts") {
    val exons1 = Iterable(Exon("e1", "t1", true, ReferenceRegion("1", 0L, 101L)),
      Exon("e2", "t1", true, ReferenceRegion("1", 200L, 401L)),
      Exon("e3", "t1", true, ReferenceRegion("1", 500L, 576L)))
    val exons2 = Iterable(Exon("e1", "t2", false, ReferenceRegion("1", 600L, 651L)),
      Exon("e2", "t2", false, ReferenceRegion("1", 200L, 401L)),
      Exon("e3", "t2", false, ReferenceRegion("1", 125L, 176L)),
      Exon("e4", "t2", false, ReferenceRegion("1", 25L, 76L)))

    val transcripts = Seq(Transcript("t1", Seq("t1"), "g1", true, exons1, Iterable(), Iterable()),
      Transcript("t2", Seq("t2"), "g1", false, exons2, Iterable(), Iterable()))
    val rdd = sc.parallelize(transcripts)

    val lengths = Quantify.extractTranscriptLengths(rdd)

    assert(lengths.size === 2)
    assert(lengths("t1") === 375L)
    assert(lengths("t2") === 350L)
  }

  def dummyTranscript(id: String): Transcript = {
    return new Transcript(id,
      Seq("test"),
      "Gene1",
      true,
      Iterable(),
      Iterable(),
      Iterable())
  }

  test("dummy transcript correctly initialized") {
    var t = dummyTranscript("t1")
    assert(t.id == "t1")
    assert(t.strand == true)
  }

  sparkTest("transcripts correctly matched with coverage") {
    var s1: Double = 1
    var s2: Double = 2
    var s3: Double = 3
    val rdd1 = sc.parallelize(Array(dummyTranscript("t1"),
      dummyTranscript("t2"),
      dummyTranscript("t3")))
    val rdd2 = sc.parallelize(Array(("t1", s1, Iterable[Long]()),
      ("t2", s2, Iterable[Long]()),
      ("t3", s3, Iterable[Long]())))
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

  sparkTest("quantify unique transcripts") {
    // generate transcripts
    val tLen = Seq(1000, 600, 400, 550, 1275, 1400)
    val (transcripts,
      names,
      kmerMap,
      classMap) = TranscriptGenerator.generateIndependentTranscripts(20,
      tLen,
      Some(1234L))

    // generate reads
    val reads = ReadGenerator(transcripts, Seq(0.2, 0.1, 0.3, 0.2, 0.1, 0.1), 10000, 75, Some(4321L))

    // run quantification
    val relativeAbundances = Quantify(
      AlignmentRecordRDD(sc.parallelize(reads), SequenceDictionary.empty, RecordGroupDictionary.empty, Seq.empty),
      sc.parallelize(kmerMap.toSeq),
      sc.parallelize(classMap.toSeq),
      sc.parallelize(names.zip(tLen).map(p => Transcript(p._1,
        Seq(p._1),
        p._1,
        true,
        Iterable(Exon(p._1 + "exon",
          p._1,
          true,
          ReferenceRegion(p._1, 0, p._2.toLong))),
        Iterable(),
        Iterable()))),
      20,
      20,
      false,
      false).collect
      .map(kv => (kv._1.id, kv._2))
      .toMap

    assert(relativeAbundances.size === 6)
    assert(fpEquals(relativeAbundances("0"), 0.2, 0.05))
    assert(fpEquals(relativeAbundances("1"), 0.1, 0.05))
    assert(fpEquals(relativeAbundances("2"), 0.3, 0.05))
    assert(fpEquals(relativeAbundances("3"), 0.2, 0.05))
    assert(fpEquals(relativeAbundances("4"), 0.1, 0.05))
    assert(fpEquals(relativeAbundances("5"), 0.1, 0.05))
  }

  sparkTest("Test of TestingTwoBitFile") {
    val region1 = ReferenceRegion("region1", 0L, 10L)
    val tbf = new TestingTwoBitFile(new ByteArrayByteAccess(new Array[Byte](1)))
    assert(tbf.extract(region1) == "CAATCCTTCG")
  }

  sparkTest("Test of Index") {
    // Takes a set of transcripts and a twobitfile and a kmer length, then returns a tuple: (kmers -> eq classes, eq class -> iterable of member kmers)

    val region1 = ReferenceRegion("region1", 0L, 10L)
    val exon1 = Exon("exon1", "transcript1", true, region1)
    val transcript1 = Transcript("transcript1", Seq("transcript1"), "gene1", true, Iterable(exon1), Iterable(), Iterable())

    val region2 = ReferenceRegion("region2", 11L, 20L)
    val exon2 = Exon("exon2", "transcript2", true, region2)
    val transcript2 = Transcript("transcript2", Seq("transcript2"), "gene1", true, Iterable(exon2), Iterable(), Iterable())

    val transcripts = sc.parallelize(Seq(transcript1, transcript2))

    val tbfile = new TestingTwoBitFile(new ByteArrayByteAccess(new Array[Byte](1)))

    // List of ( ... (kmer, class id) ... ) AND (... (id, list of kmers) ...) 
    var (kmersToEq, eqToKmers) = Index.apply(tbfile, transcripts, 5)

    val kToEq = kmersToEq.collect()
    val eqToK = eqToKmers.collect()

    // Tracking 3 particular kmers: "CAATC", "GTGCA", "CTTCG"
    // Should be at least 2 eq classes, such that one contains both "CAATC" and "CTTCG", while the ther contains "GTGCA"
    val CAATCtoClassArray = kToEq.filter(_._1 == "CAATC")
    assert(CAATCtoClassArray.length == 1) // There should only be one instance of each kmer
    val CAATCtoClass = CAATCtoClassArray(0)

    val GTGCAtoClassArray = kToEq.filter(_._1 == "GTGCA")
    assert(GTGCAtoClassArray.length == 1)
    val GTGCAtoClass = GTGCAtoClassArray(0)

    val CTTCGtoClassArray = kToEq.filter(_._1 == "CTTCG")
    assert(CTTCGtoClassArray.length == 1)
    val CTTCGtoClass = CTTCGtoClassArray(0)

    var class1 = CAATCtoClass._2
    var class2 = GTGCAtoClass._2
    var class3 = CTTCGtoClass._2
    assert(class1 != class2)
    assert(class1 == class3)

    var class1Kmers = eqToK.filter(_._1 == class1)(0)
    assert(class1Kmers._2.toList.contains("CAATC"))
    assert(class1Kmers._2.toList.contains("CTTCG"))
    assert(!class1Kmers._2.toList.contains("GTGCA"))

    var class2Kmers = eqToK.filter(_._1 == class2)(0)
    assert(class2Kmers._2.toList.contains("GTGCA"))
    assert(!class2Kmers._2.toList.contains("CAATC"))
    assert(!class2Kmers._2.toList.contains("CTTCG"))
  }

  sparkTest("quantify a small set of more realistic but unbiased transcripts") {
    // generate transcripts
    val classSize = Seq(1000, 500, 700, 400, 400, 200, 100)
    val classMultiplicity = Seq(1, 1, 1, 1, 2, 2, 3)
    val classMembership = Seq(Set(0),
      Set(1, 2),
      Set(1, 3),
      Set(1, 4),
      Set(2, 5),
      Set(2, 6),
      Set(3, 6),
      Set(6))
    val (transcripts,
      names,
      kmerMap,
      classMap) = TranscriptGenerator.generateTranscripts(20,
      classSize,
      classMultiplicity,
      classMembership,
      Some(1000L))
    val tLen = transcripts.map(_.length)

    // generate reads
    val abundances = Seq(0.05, 0.1, 0.25, 0.1, 0.05, 0.025, 0.025, 0.4)
    val reads = ReadGenerator(transcripts,
      abundances,
      50000,
      75,
      Some(5000L))

    // run quantification
    val relativeAbundances = Quantify(
      AlignmentRecordRDD(sc.parallelize(reads), SequenceDictionary.empty, RecordGroupDictionary.empty, Seq.empty),
      sc.parallelize(kmerMap.toSeq),
      sc.parallelize(classMap.toSeq),
      sc.parallelize(names.zip(tLen).map(p => Transcript(p._1,
        Seq(p._1),
        p._1,
        true,
        Iterable(Exon(p._1 + "exon",
          p._1,
          true,
          ReferenceRegion(p._1, 0, p._2.toLong))),
        Iterable(),
        Iterable()))),
      20,
      50,
      false,
      false).collect
      .map(kv => (kv._1.id, kv._2))
      .toMap

    assert(relativeAbundances.size === 8)
    assert(fpEquals(relativeAbundances("0"), 0.05, 0.01))
    assert(fpEquals(relativeAbundances("1"), 0.1, 0.05))
    assert(fpEquals(relativeAbundances("2"), 0.25, 0.05))
    assert(fpEquals(relativeAbundances("3"), 0.1, 0.05))
    assert(fpEquals(relativeAbundances("4"), 0.05, 0.025))
    assert(fpEquals(relativeAbundances("5"), 0.025, 0.0125))
    assert(fpEquals(relativeAbundances("6"), 0.025, 0.0125))
    assert(fpEquals(relativeAbundances("7"), 0.4, 0.05))
  }

  sparkTest("quantify unique transcripts without bias removal") {
    // generate transcripts
    val tLen = Seq(1000, 600, 400, 550, 1275, 1400)
    val (transcripts,
      names,
      kmerMap,
      classMap) = TranscriptGenerator.generateIndependentTranscripts(20,
      tLen,
      Some(1234L))

    // generate reads
    val reads = ReadGenerator(transcripts, Seq(0.2, 0.1, 0.3, 0.2, 0.1, 0.1), 10000, 75, Some(4321L))

    // run quantification
    val relativeAbundances = Quantify(
      AlignmentRecordRDD(sc.parallelize(reads), SequenceDictionary.empty, RecordGroupDictionary.empty, Seq.empty),
      sc.parallelize(kmerMap.toSeq),
      sc.parallelize(classMap.toSeq),
      sc.parallelize(names.zip(tLen).map(p => Transcript(p._1,
        Seq(p._1),
        p._1,
        true,
        Iterable(Exon(p._1 + "exon",
          p._1,
          true,
          ReferenceRegion(p._1, 0, p._2.toLong))),
        Iterable(),
        Iterable()))),
      20,
      20,
      false,
      false).collect
      .map(kv => (kv._1.id, kv._2))
      .toMap

    assert(relativeAbundances.size === 6)
    assert(fpEquals(relativeAbundances("0"), 0.2, 0.05))
    assert(fpEquals(relativeAbundances("1"), 0.1, 0.05))
    assert(fpEquals(relativeAbundances("2"), 0.3, 0.05))
    assert(fpEquals(relativeAbundances("3"), 0.2, 0.05))
    assert(fpEquals(relativeAbundances("4"), 0.1, 0.05))
    assert(fpEquals(relativeAbundances("5"), 0.1, 0.05))
  }

  sparkTest("quantify unique transcripts where length bias is so strong that all variation in transcript abundance is due to length") {
    // generate transcripts
    val tLen = Seq(1000, 600, 400, 550, 1275, 1400)
    val (transcripts,
      names,
      kmerMap,
      classMap) = TranscriptGenerator.generateIndependentTranscripts(20,
      tLen,
      Some(1234L))
    // generate reads
    val totLen = tLen.sum.toDouble
    val reads = ReadGenerator(transcripts, tLen.map(x => x / totLen), 10000, 75, Some(4321L))

    // run quantification
    val relativeAbundances = Quantify(
      AlignmentRecordRDD(sc.parallelize(reads), SequenceDictionary.empty, RecordGroupDictionary.empty, Seq.empty),
      sc.parallelize(kmerMap.toSeq),
      sc.parallelize(classMap.toSeq),
      sc.parallelize(names.zip(tLen).map(p => Transcript(p._1,
        Seq(p._1),
        p._1,
        true,
        Iterable(Exon(p._1 + "exon",
          p._1,
          true,
          ReferenceRegion(p._1, 0, p._2.toLong))),
        Iterable(),
        Iterable()))),
      20,
      20,
      true,
      true).collect
      .map(kv => (kv._1.id, kv._2))
      .toMap

    assert(relativeAbundances.size === 6)
    assert(fpEquals(relativeAbundances("0"), 1.0 / 6, 0.05))
    assert(fpEquals(relativeAbundances("1"), 1.0 / 6, 0.05))
    assert(fpEquals(relativeAbundances("2"), 1.0 / 6, 0.05))
    assert(fpEquals(relativeAbundances("3"), 1.0 / 6, 0.05))
    assert(fpEquals(relativeAbundances("4"), 1.0 / 6, 0.05))
    assert(fpEquals(relativeAbundances("5"), 1.0 / 6, 0.05))
  }

  sparkTest("quantify unique transcripts with a weaker length bias") {
    // generate transcripts
    val tLen = Seq(1000, 600, 400, 550, 1275, 1400) // average length is 870.83
    val (transcripts,
      names,
      kmerMap,
      classMap) = TranscriptGenerator.generateIndependentTranscripts(20,
      tLen,
      Some(1234L))

    // generate reads
    val reads = ReadGenerator(transcripts, Seq(0.2, 0.1, 0.05, 0.2, 0.05, 0.4), 10000, 75, Some(4321L))

    // run quantification
    val relativeAbundances = Quantify(
      AlignmentRecordRDD(sc.parallelize(reads), SequenceDictionary.empty, RecordGroupDictionary.empty, Seq.empty),
      sc.parallelize(kmerMap.toSeq),
      sc.parallelize(classMap.toSeq),
      sc.parallelize(names.zip(tLen).map(p => Transcript(p._1,
        Seq(p._1),
        p._1,
        true,
        Iterable(Exon(p._1 + "exon",
          p._1,
          true,
          ReferenceRegion(p._1, 0, p._2.toLong))),
        Iterable(),
        Iterable()))),
      20,
      20,
      true,
      true).collect
      .map(kv => (kv._1.id, kv._2))
      .toMap

    assert(relativeAbundances.size === 6)

    // Transcript 2 is the shortest at 440 bp, and is one of the least abundant
    // Part of this low abundance is bias due to short length, so calibration
    // should increase its abundance.
    assert(relativeAbundances("2") > 0.05)

    // Transcript 5 is the longest at 1400 bp. It is the most abundant.
    // Part of this high abudnace is bias due to long length.
    // Calibration should therefore decrease its abundance.
    assert(relativeAbundances("5") < 0.4)
  }

}
