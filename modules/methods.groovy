
/*
 * nourdine.bah@crick.ac.uk
 */

import java.nio.file.Paths

@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv

/////////////////////
def absPath(path) {//
/////////////////////

	def f = new File(path)

	if ( ! f.isAbsolute() ) {
		return Paths.get(workflow.projectDir.toString(), path).toString()
	} else {
		return path
	}
}

/////////////////////////////////
def addValue(map, key, value) {//
/////////////////////////////////
	def new_map = map.clone()
	new_map.put(key, value)
	return new_map
}

///////////////////////////////////
def parseSeries(metadata, path) {//
///////////////////////////////////
	def csv = parseCsv( new File(path).text )
	def channels = []
	for ( row in csv ) {
		channels.add( metadata.clone() << row.toMap() )
	}
	return channels
}

///////////////////////////
def dropKeys(map, keys) {//
///////////////////////////
	// because Map.dropWhile doesn't with the current Java version of Nextflow,
	// apparently requires Java9
	def new_map = [:]
	map.each{ k, v ->
		if ( ! keys.contains(k) ) {
			new_map[k] = v
		}
	}
	return new_map
}


