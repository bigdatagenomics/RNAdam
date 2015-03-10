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

import org.apache.spark.SparkContext._
import org.apache.spark.rdd.RDD
import org.apache.spark.broadcast.Broadcast
import org.apache.spark.mllib.linalg.Vectors
import org.apache.spark.mllib.regression.{ LabeledPoint, LinearRegressionModel, LinearRegressionWithSGD }
import scala.math.{ exp, log }
import org.jblas.{ DoubleMatrix, Solve }

object Tare extends Serializable {

  /**
   * Converts a base into an integer number. Only works for A/C/G/T, both
   * lower and upper case. Will throw an exception for other bases, user is
   * expected to check before passing to function.
   *
   * @param base The base to convert into an integer.
   * @return Returns a mapping from a base to an integer.
   */
  private def idx(base: Char): Int = base match {
    case 'A' | 'a' => 0
    case 'C' | 'c' => 1
    case 'G' | 'g' => 2
    case 'T' | 't' => 3
  }

  /**
   * @param base Base to examine.
   * @return Returns true if the base is a "valid" base as defined by the idx function.
   *
   * @see idx
   */
  private def validBase(base: Char): Boolean = base match {
    case 'A' | 'a' | 'C' | 'c' | 'G' | 'g' | 'T' | 't' => true
    case _ => false
  }

  /**
   * Hashes a dinucleotide sequence context into an integer.
   *
   * @param context Sequence context to hash.
   * @return Mapping from sequence context into an int; scale is 0-15.
   */
  private def dinucToIdx(context: String): Int = {
    4 * idx(context(0)) + idx(context(1))
  }

  /**
   * @param context Dinucleotide sequence context to evaluate.
   * @return Returns true if the context is valid. Validity is defined by the
   *         idx function.
   *
   * @see idx
   */
  private def isValidContext(context: String): Boolean = {
    context.length == 2 &&
      validBase(context(0)) &&
      validBase(context(1))
  }

  /**
   * Featurizes a k-mer. Cuts the k-mer into dinucleotide sequence contexts, filters out invalid
   * sequences, and then builds a feature vector that maps the log of the multiplicity to the
   * fraction of each context occurrence.
   *
   * @param kmer K-mer string to featurize.
   * @param multiplicity Observed k-mer count.
   * @return Vector mapping k-mer to log of multiplicity and sequence context counts.
   */
  private[quantification] def kmerToDinucFeature(kmer: String, multiplicity: Long): LabeledPoint = {
    // slice kmers into contexts
    val contexts = kmer.sliding(2).filter(isValidContext).toIterable
    assert(contexts.size > 0, "k-mer: " + kmer + " does not contain any valid contexts.")
    val ctxReciprocal = 1.0 / contexts.size.toDouble

    // populate initial context array
    val contextCounts = Array.fill[Double](16) { 0.0 }

    // update counts
    contexts.foreach(c => contextCounts(dinucToIdx(c)) += ctxReciprocal)

    new LabeledPoint(log(multiplicity.toDouble), Vectors.dense(contextCounts))
  }

  /**
   * Translates k-mers and counts into features and then trains a linear regression model.
   * After the model is trained, the k-mer count is recalibrated.
   *
   * @param kmers An RDD of k-mer strings and counts.
   * @return Returns a calibrated RDD of k-mers and counts.
   */
  def calibrateKmers(kmers: RDD[(String, Long)]): RDD[(String, Long)] = {
    // featurize kmers
    val countAccumulator = kmers.context.accumulator(0L)
    val multiplicityAccumulator = kmers.context.accumulator(0L)
    val kmersAndFeatures = kmers.map(kv => {
      countAccumulator += 1
      multiplicityAccumulator += kv._2
      (kv._1, kmerToDinucFeature(kv._1, kv._2))
    }).cache()

    // train linear model
    val model = new LinearRegressionWithSGD().setIntercept(true).run(kmersAndFeatures.map(kv => kv._2))
    val bcastModel = kmersAndFeatures.context.broadcast(model)

    // compute sample mean
    val mean = log(multiplicityAccumulator.value.toDouble / countAccumulator.value.toDouble)

    // transform back 
    val calibratedKmers = kmersAndFeatures.map(kv => {
      (kv._1, exp(mean + (kv._2.label - bcastModel.value.predict(kv._2.features))).toLong)
    })

    // unpersist feature rdd
    kmersAndFeatures.unpersist()

    calibratedKmers
  }

  /**
   * Trains a linear regression model of the log of the transcript abundances
   * versus the log of the transcript lengths.
   * This is used to recalibrate the transcript abundances.
   * Warning: This procedure assumes that the abundances are all > 0
   * The log is undefined for numbers <= 0
   *
   * @param muHat An RDD containing the transcript abundances.
   * @param tLen A mapping from transcript names to their lengths.
   * @param samplingRate The fraction of the RDD that should be picked a sample when performing linear regression.
   * @return The recalibrated transcript abundances.
   */
  def calibrateTxLenBias(muHat: RDD[(String, Double, Iterable[Long])],
                         tLen: scala.collection.Map[String, Long],
                         samplingRate: Double = 1.0): RDD[(String, Double, Iterable[Long])] = {
    // A local (not distributed) copy of muHat without the equivalence
    // class ids. Removing these ids saves space since they are
    // not used here.
    val local: Array[(String, Double)] = muHat.map(mh => (mh._1, mh._2)).sample(false, samplingRate).collect()

    // Calculate the log of the average (mean) abundance.
    // Since all the abundances sum to 1.0 by definition,
    // this is simply 1.0 divided by the number of transcripts.
    // Use the fact that log(1.0 / x) = -log(x)
    val mean: Double = -log(local.length.toDouble)

    // Do linear regression by constructing and solving matrix equation.
    // We did not use Spark MLLib because it did not converge the correct
    // values. The problem appears to be with using one-dimensional feature
    // vectors on small sample sizes. So we used JBLAS Linear Algebra for Java.
    val y: DoubleMatrix = new DoubleMatrix(local.map(mh => log(mh._2)))
    val slopev: DoubleMatrix = new DoubleMatrix(local.map(mh => log(tLen(mh._1).toDouble)))
    val intv: DoubleMatrix = DoubleMatrix.ones(local.length)
    val x: DoubleMatrix = DoubleMatrix.concatHorizontally(slopev, intv)
    val xt: DoubleMatrix = x.transpose()
    val xtx: DoubleMatrix = xt.mmul(x)
    val xty: DoubleMatrix = xt.mmul(y)
    val soln: DoubleMatrix = Solve.solve(xtx, xty)
    val slope: Double = soln.get(0, 0)
    val intercept: Double = soln.get(1, 0)

    // Use the slope and the intercept found through regression
    // to calibrate transcript abundances to correct for variation
    // in  transcript lengths. These calibrated abundances may not
    // add up to 1 because regression was done on a logarithmic scale.
    // muHat is used here instead of local because in the future, local
    // might become a sample of muHat instead of a full copy.
    // Furthermore, doing too much locally may create a larger than
    // necessary bottleneck in the computation.
    val calMuHat: RDD[(String, Double, Iterable[Long])] = muHat.map(mh => (mh._1, exp(mean + (((slope * mh._2) + intercept) - mh._2)), mh._3))

    val totalAbundance: Double = calMuHat.map(mh => mh._2).reduce((x: Double, y: Double) => x + y)

    // Return the calibrated abundances.
    calMuHat.map((mh: (String, Double, Iterable[Long])) => (mh._1, mh._2 / totalAbundance, mh._3))
  }
}
