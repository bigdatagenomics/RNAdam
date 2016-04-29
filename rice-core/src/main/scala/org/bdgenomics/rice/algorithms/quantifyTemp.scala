/*// TEMP QUANTIFY:

// SHOULDNT LIKELIHOOD BE A FLOAT ?
def likelihoodEstimator(vals: Iterable[String, Long]): Map[String, Double] = {
	// Total number of occurrences of the kmer
	val total = vals.reduce( (v1, v2) => (0, v1._2 + v2._2) )._2

	// Likelihood is the occurrences in a certain transcript as a proportion of the total 
	val likelihood = vals.map( v => (v._1, {v._2.toDouble / total.toDouble})) 

}

def prep(reads: RDD[AlignmentRecord],
 		 graph: Graph[V, E]): RDD[(Long, Map[String, Double])] = {

	// Prep Vertices:
	vertices = graph.vertices 													  // RDD[ kmerHash, ColoredKmerVertex ]
 					.flatmap(v => [ (v._1, (t._1, 1L)) for t in v._2.terminals ]) // RDD[ kmerHash, (color, 1) ]
 					.reduceByKey( (c1, c2) => (c1._1, c1._2 + c2._2) ) 			  // RDD[ kmerHash, (color, num occurrences) ]   ** REPLACE WITH SINGLE MAP
 					.groupByKey() 												  // RDD[ kmerHash, Iterable[(color, num occurrences)] ]
 					.map( v => likelihoodEstimator(v) ) 						  // RDD[ kmerHash, Map[color -> likelihood] ] ** MOVE TO AFTER READ DATA IS ADDED

 	// Prep Edges:
 	preppedReads = reads.zipWithUniqueId()																// RDD[ read, readID ]
 						.map(r => (r._1.sequence, r._2))												// RDD[ sequence, readID ]
 						.flatmap( r => [ (imer.longHash, r._2) for imer in Intmer.fromSequence(r._1) ]) // RDD[ kmerHash, readID]

 	// Join:
 	preppedReads.join(vertices) // RDD[ kmerHash, ( readID, Map[color -> likelihood] ) ]
 				.map(v => v._2) // RDD[ readId, Map[color -> likelihood] ]
}


  */ 