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
import org.apache.spark.SparkContext._
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.formats.avro.NucleotideContigFragment
import org.bdgenomics.rice.Timers._
import org.apache.spark.graphx.Graph
import net.fnothaft.ananas.debruijn.{ColoredDeBruijnGraph, ColoredKmerVertex}
import net.fnothaft.ananas.models.ContigFragment

object Index extends Serializable with Logging {

  /**
   * Computes an index, given a set of Nucleotide Contig Fragments. An index provides a de bruijn graph of kmers
   *
   * @param contigFragments An RDD containing contigFragments.
   * @return Returns a Graph representing a colored De Bruijn graph of kmers
   */
  def apply(contigFragments: RDD[ContigFragment]): Graph[ColoredKmerVertex, Unit] = {

    ColoredDeBruijnGraph.buildFromFragments(contigFragments)

  }

}