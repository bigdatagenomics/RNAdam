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
package org.bdgenomics.rice.models

import net.fnothaft.ananas.models.CanonicalKmer
import java.io.{ ObjectInputStream, FileInputStream }

object KmerIndex extends Serializable {

  def apply(filename: String): KmerIndex = {
    val in = new ObjectInputStream(new FileInputStream(filename))
    val mapping = in.readObject().asInstanceOf[Map[(Long, Boolean), Map[String, Long]]]
    in.close()
    IndexMap(16, mapping)
  }
}

trait KmerIndex extends Serializable {

  def getKmerLength: Int

  def getTranscripts(kmer: CanonicalKmer): Map[String, Long]
}

case class IndexMap(kmerLength: Int,
                    kmersToCounts: Map[(Long, Boolean), Map[String, Long]]) extends KmerIndex {

  def getKmerLength: Int = kmerLength

  def getTranscripts(kmer: CanonicalKmer): Map[String, Long] = kmersToCounts((kmer.longHash, kmer.isOriginal))

}
