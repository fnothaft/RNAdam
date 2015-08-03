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
import org.bdgenomics.formats.avro.{ Contig, NucleotideContigFragment }
import org.bdgenomics.adam.models.{ Exon, ReferenceRegion, Transcript, CDS, UTR }
import org.bdgenomics.rice.utils.riceFunSuite
import scala.collection.Map
import scala.collection.immutable.HashMap
import scala.math.abs
import org.bdgenomics.utils.io.{ ByteAccess, ByteArrayByteAccess }
import net.fnothaft.ananas.models.{CanonicalKmer, IntMer, ContigFragment}
import net.fnothaft.ananas.avro.{Kmer, Backing}
import net.fnothaft.ananas.debruijn.ColoredDeBruijnGraph

class QuantifySuite extends riceFunSuite {

  // Blatant copy of function from ananas/models/ContigFragment to get around protected status
  def def buildFromNCF(fragment: NucleotideContigFragment): ContigFragment = {
    // is this the last fragment in a contig?
    val isLast = fragment.getFragmentNumber == (fragment.getNumberOfFragmentsInContig - 1)

    val sequence = IntMer.fromSequence(fragment.getFragmentSequence)
      .map(_.asInstanceOf[CanonicalKmer])

    new ContigFragment(fragment.getContig.getContigName,
                       sequence,
                       isLast,
                       Option(fragment.getFragmentStartPosition).fold(0)(_.toInt))
  }

  def createTestIndex(sequence: String = "ACACTGTGGGTACACTACGAGA") : (Map[Long, Map[String, Long]], Map[String, Transcript]) = {
    val ncf = NucleotideContigFragment.newBuilder()
      .setContig(Contig.newBuilder()
      .setContigName("ctg")
      .build())
      .setFragmentNumber(0)
      .setNumberOfFragmentsInContig(1)
      .setFragmentSequence(sequence)
      .build()

    val frag = sc.parallelize( Seq( buildFromNCF(ncf) ) ) 

    val tx = Seq( Transcript("one", Seq("one"), "gene1", true, Iterable[Exon](), Iterable[CDS](), Iterable[UTR]()) ,
                  Transcript("two", Seq("two"), "gene1", true, Iterable[Exon](), Iterable[CDS](), Iterable[UTR]()) )
    val transcripts = sc.parallelize(tx)

    Index(frag, transcripts)
  }

  sparkTest("Simple Test of Index") {
    val testSeq = "ACACTGTGGGTACACTACGAGA"
    val (imap, tmap) = createTestIndex(testSeq)

    // Test kmer mapping
    val imers = testSeq.sliding(16).map(s => IntMer(s))
    assert( imap.size == 7 ) // 7 kmers of length 16
    assert( imers.forall(i => imap(i.longHash)("ctg") == 1) )

    // Test transcript mapping
    assert( tmap("one").id == "one")
    assert( tmap("two").id == "two")
  }
}
