import groovy.json.JsonSlurper
import groovy.transform.CompileDynamic
import groovy.transform.Immutable
import java.nio.file.Path
import nextflow.Nextflow
import nextflow.io.ValueObject
import nextflow.util.KryoHelper

@CompileDynamic
@ValueObject
@Immutable(copyWith = true)
class References extends Nextflow {
    static {
        KryoHelper.register(References)
    }

    static References load(Path json, String annotationName) {
        JsonSlurper jsonSlurper = new JsonSlurper()

        Map<String, Path> references = jsonSlurper
            .parse(json)
            .collectEntries { name, path -> [name, file(path)] }
        
        return new References(references, annotationName)
    }

    /** Mapping of reference name to path */
    Map<String, Path> references
    /** The name of the reference to use for annotation */
    String annotationName

    Path dir(String name) {
        if (!references.containsKey(name)) {
            throw new RuntimeException("Invalid reference name ${name}")
        }
        return references[name]
    }

    Path getAnnotationDir() {
        return dir(annotationName)
    }

    String getAnnotationFastaName() {
        return "${annotationName}.fasta.gz"
    }
}