//
// This file holds several Groovy functions that could be useful for any Nextflow pipeline
//

import org.yaml.snakeyaml.Yaml
import java.nio.file.Paths

@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv

class Utils {
    //
    // When running with -profile conda, warn if channels have not been set-up appropriately
    //
    public static void checkCondaChannels(log) {
        Yaml parser = new Yaml()
        def channels = []
        try {
            def config = parser.load("conda config --show channels".execute().text)
            channels = config.channels
        } catch(NullPointerException | IOException e) {
            log.warn "Could not verify conda channel configuration."
            return
        }

        // Check that all channels are present
        def required_channels = ['conda-forge', 'bioconda', 'defaults']
        def conda_check_failed = !required_channels.every { ch -> ch in channels }

        // Check that they are in the right order
        conda_check_failed |= !(channels.indexOf('conda-forge') < channels.indexOf('bioconda'))
        conda_check_failed |= !(channels.indexOf('bioconda') < channels.indexOf('defaults'))

        if (conda_check_failed) {
            log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  There is a problem with your Conda configuration!\n\n" +
                "  You will need to set-up the conda-forge and bioconda channels correctly.\n" +
                "  Please refer to https://bioconda.github.io/user/install.html#set-up-channels\n" +
                "  NB: The order of the channels matters!\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        }
    }

    //
    // NOURDINE CUSTOM FUNCTIONS
    //

    public static void absPath(path) {
        def f = new File(path)

        if ( ! f.isAbsolute() ) {
            return Paths.get(workflow.projectDir.toString(), path).toString()
        } else {
            return path
        }
    }

    public static void addValue(map, key, value) {
        def new_map = map.clone()
        new_map.put(key, value)
        return new_map
    }

    public static void parseSeries(metadata, path) {
        def csv = parseCsv( new File(path).text )
        def channels = []
        for ( row in csv ) {
            channels.add( metadata.clone() << row.toMap() )
        }
        return channels
    }

    public static void dropKeys(map, keys) {
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

    public static void removeKeys(map, keys) {
        def new_map = [:]

        map.each{
            if ( ! keys.contains(it.key) )
            {
                new_map.put(it.key, it.value)
            }
        }

        return new_map
    }

    public static void getMinLength(structure) {
        return structure.split("[A-Z]").collect{it as int }.sum()
    }

    public static void getPuckName(puck) {
        def f = new File(puck)
        return f.getName().toString().replaceAll('\\.csv$', '') // single quotes!
    }
}
