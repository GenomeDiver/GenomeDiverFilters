package recommend 

import java.io.File
import scopt.OParser
import zio.logging.console
import scala.jdk.CollectionConverters._
import collection.convert.ImplicitConversions._
import scala.jdk.OptionConverters._

import org.monarchinitiative.phenol.ontology.data.{Term, TermId}
import org.monarchinitiative.phenol.annotations.io.hpo.DiseaseDatabase._
import org.monarchinitiative.phenol.annotations.io.hpo._
import org.monarchinitiative.phenol.annotations.formats.hpo.{HpoAssociationData, HpoDisease}
import org.monarchinitiative.phenol.io.OntologyLoader
import org.monarchinitiative.lirical.io.LiricalDataResolver
import org.monarchinitiative.lirical.core.service.PhenotypeService
import org.monarchinitiative.lirical.core.likelihoodratio.{InducedDiseaseGraph, PhenotypeLikelihoodRatio}

import java.io.{File}
import java.nio.file.{Paths, Path}
import java.io.PrintWriter

// import org.virtuslab.iskra.api._
import org.apache.spark.sql.SparkSession
import org.apache.spark.sql.expressions.Window
import org.apache.spark.sql.functions._
import org.apache.spark.sql.types.{IntegerType,StringType,StructType,StructField}
import org.apache.spark.sql.catalyst.ScalaReflection


object PhenotypeRecommendation {

  case class Config(
    hpo:Seq[String] = Seq(), 
    num:Int = 20, 
    sex:String = "UNKNOWN",
    lirical_data:File = new File("."),
    genes_tsv:File = new File("."),
    ontology:File = new File("."),
    annotation:File = new File("."), 
    output:File = new File(".")
  )
  
  val gd_version = "1.0.0" 
  
  val builder = OParser.builder[Config]

  val arg_parser = {
    import builder._

    OParser.sequence(
      programName("Phenotype Recommendation"),
      head("Genome Diver", gd_version),
      opt[Int]('e', "num")
        .optional()
        .valueName("<num>")
        .action((x, c) => c.copy(num = x))
        .text("num is optional. It is the number of recommended phenotypes (default: 20)"),
      opt[Seq[String]]('h', "hpo")
        .required().valueName("<hpo1>, <hpo2>")
        .action((x, c) => c.copy(hpo = x))
        .text("hpo is a required property. It is an input of Exomiser"),
      opt[String]('s', "sex")
        .optional()
        .valueName("<sex>")
        .action((x, c) => c.copy(sex = x))
        .text("sex is optional. [MALE, FEMALE, UNKNOWN] Default is UNKNOWN"),
      opt[File]('g', "genes_tsv")
        .required()
        .valueName("<file>")
        .action((x, c) => c.copy(genes_tsv = x))
        .text("genes_tsv is a required file property. It is an output of Exomiser."),
      opt[File]('l', "lirical_data")
        .required()
        .valueName("<directory>")
        .action((x, c) => c.copy(lirical_data = x))
        .text("lirical_data is a required directory. It is an input of the phenotype recommendation."), 
      opt[File]('o', "output")
        .required()
        .valueName("<file>")
        .action((x, c) => c.copy(output = x))
        .text("output is a required file property. This is the output file (txt format).")
    )
  }
  
  def phenotype_recommendation(config: Config) : Unit = {
    implicit val spark: SparkSession = SparkSession
      .builder().master("local")
      .appName("gd").getOrCreate()
    
    import spark.implicits._
    spark.sparkContext.setLogLevel("ERROR")
    
    // Nathan proposes liklihood cutoff. 
    // NUM_RECOMMEND .. 
    // Causative Gene may not be in top 4, try 10 
    
    /* hardcoded mess 
    val NUM_RECOMMEND = 20
    val LIR_PATH = Paths.get("/home/kshi/gd-data/lirical-cli-2.0.0-RC2-json/data")
    val ONT_PATH = Paths.get("/home/kshi/gd-data/lirical-cli-2.0.0-RC2-json/data/hp.json")
    val ANN_PATH = Paths.get("/home/kshi/gd-data/lirical-cli-2.0.0-RC2-json/data/phenotype.hpoa")
    val INPUT_HPO = List("HP:0001156","HP:0001363","HP:0010055","HP:0011304")
    val GENES_TSV = Paths.get("./patient/patient-1/analysis-1/rare.genes.tsv")
    val OUTPUT_HPO = Paths.get("./patient/patient-1/analysis-1/recommend.txt")
    */
    
    val NUM_RECOMMEND = config.num
    val LIR_PATH = config.lirical_data.toPath
    val INPUT_HPO = config.hpo
    val SEX = config.sex
    val GENES_TSV = config.genes_tsv.toPath
    val OUTPUT_HPO = config.output.toPath
    
    val liricalDataResolver = LiricalDataResolver.of(LIR_PATH)
    val custom_options = HpoDiseaseLoaderOptions.of(
      Set(DiseaseDatabase.OMIM, DiseaseDatabase.DECIPHER).asJava, true,
      HpoDiseaseLoaderOptions.DEFAULT_COHORT_SIZE)

    val ontology = OntologyLoader.loadOntology(liricalDataResolver.hpoJson.toFile)
    val loader = HpoDiseaseLoaders.v2(ontology, HpoDiseaseLoaderOptions.defaultOptions());
    val diseases = loader.load(liricalDataResolver.phenotypeAnnotations)
    val associationData = HpoAssociationData.builder(ontology)
      .hgncCompleteSetArchive(liricalDataResolver.hgncCompleteSet())
      .mim2GeneMedgen(liricalDataResolver.mim2geneMedgen())
      .hpoDiseases(diseases)
      .build()

    val phenotypeService = PhenotypeService.of(ontology, diseases, associationData)
    val disease_map = diseases.diseaseById()
    val term_map = ontology.getTermMap
    val pheno_lr = new PhenotypeLikelihoodRatio(ontology, diseases)
     
    /* 
     * Phenotype Recommendation 
     * =============================
     * + Gather disease (ORPHA | OMIM) from top scoring genes (via Exomiser)
     * + Use HPO Annotations to find affiliated HPO phenotypes 
     * - Pare down Phenotypes against input Phenotypes
     * - Further pare down Phenotypes against input->ancestors 
     *
     * */
   
    case class Recommendation(
      hpo_id: String,
      hpo_name: String,
      hpo_definition: String,
      disease: String, 
      lr_scores: String,
      cv: Float
    )

    val genes = spark.read
      .option("delimiter", "\t")
      .option("header", true)
      .csv(GENES_TSV.toAbsolutePath.toString)

    val schema = new StructType()
      .add("DatabaseId",StringType,true)
      .add("DB_Name",StringType,false)
      .add("Qualifier",StringType,false)
      .add("HPO_ID",StringType,false)
      .add("DB_Reference",StringType,false)
      .add("Evidence",StringType,false)
      .add("Onset",StringType,false)
      .add("Frequency",StringType,false)
      .add("Sex",StringType,false)
      .add("Modifier",StringType,false)
      .add("Aspect",StringType,false)
      .add("BiocurationBy",StringType,false)

    val hpoa = spark.read
      .option("delimiter", "\t")
      .option("header", "false")
      .option("comment", "#")
      .schema(schema)
      .csv(LIR_PATH.resolve("phenotype.hpoa").toAbsolutePath.toString)

    hpoa.show()
    
    val top_genes = genes.select("GENE_SYMBOL", "MOI", "EXOMISER_GENE_COMBINED_SCORE", 
      "EXOMISER_GENE_PHENO_SCORE", "EXOMISER_GENE_VARIANT_SCORE", "HUMAN_PHENO_EVIDENCE", 
      "HUMAN_PPI_EVIDENCE").limit(10)
    
    println("\nExomiser Output:") 
    top_genes.show()
    
    // Diseases are taken from HUMAN_PHENO_EVIDENCE column of top genes  
    val disease_guess = top_genes.select("HUMAN_PHENO_EVIDENCE").na.drop
      .collect.flatMap(
        row => """((?:ORPHA|OMIM):\d+)""".r.findAllIn(row.getString(0))
      ).toList
    
    println("Disease Guesses:")
    println(disease_guess)

    // All phenotypes associated with diseases 
    val possible_hpo:Set[String] = disease_guess.foldLeft(Set.empty[String])(
      (acc, db_id) => {
        val annotation = disease_map.get(TermId.of(db_id)).annotationStream.iterator.asScala.toList
        val filtered_by_sex = annotation.map(_.id.toString).toSet.filter((hpo_id) => {
            val counts = hpoa.select("*")
              .where(col("DatabaseId") === db_id)
              .where(col("HPO_ID") === hpo_id)
              .where(col("Sex").isin(List(SEX):_*) || col("Sex").isNull)
              .count()
          
            counts > 0
        })
        acc.union(filtered_by_sex)
      }).diff(
        // - hpo already specified 
        // - the ancestors of specified
        INPUT_HPO.foldLeft(INPUT_HPO.toSet)( 
        (acc, hpo_id) => {
          acc ++ ontology.getAncestorTermIds(TermId.of(hpo_id))
            .asScala.toSeq.map(_.toString).toSet
        }
    ))

    val idg_cache = disease_guess.foldLeft(Map[String, InducedDiseaseGraph]()){
      (m, d) => m + (d -> InducedDiseaseGraph.create(
        disease_map.get(TermId.of(d)),
      phenotypeService.hpo()))
    }

    // TODO: investigate LR,. should be -3 .. +3
    // ... manually run through lirical-cli 
    val data = possible_hpo.toSeq.map(hpo_id => {
      val scores = disease_guess.map(disease_id =>
        pheno_lr.lrForObservedTerm(
          TermId.of(hpo_id),
          idg_cache(disease_id)
        ).lr())
      
      val best_guess = disease_guess.zip(scores).sortBy(_._2)(Ordering.Double.IeeeOrdering).last._1
      val vector = breeze.linalg.Vector[Double](scores.toArray)
      val mean = breeze.stats.mean(vector)
      val vari = breeze.stats.variance(vector)
      val term = term_map.get(TermId.of(hpo_id))
      
      // Index of Dispersion - Variance / Mean 
      val cv = (vari / mean)
      
      (
        hpo_id,
        term.getName,
        term.getDefinition,
        best_guess,
        scores.map(v => f"$v%1.3f").mkString(", "),
        cv.toFloat
      )
      }
    ).sortBy(r => (r._6))(Ordering.Float.IeeeOrdering.reverse).toDF(
      "hpo_id", "hpo_name", "hpo_definition", 
      "disease", "lr_scores", "cv")

    println("\nLR Score all phenotypes and calculate Relative Variance")
    data.show()
    
    val byDisease = Window.partitionBy("disease").orderBy(col("cv").desc)
    val disease = data.select("hpo_id", "disease", "cv")
      .withColumn("rank", rank().over(byDisease))
      .where(col("disease").isin(disease_guess:_*))
      .where(col("rank").leq((NUM_RECOMMEND / disease_guess.length).max(1)))
    
    println("\nDisease w/ relative variance rank") 
    disease.show()
 
    val output = disease.select("hpo_id").map(r => r.getString(0)).collect.toList
    val pw = new PrintWriter(new File(OUTPUT_HPO.toString))

    pw.write(output.mkString("\n"))
    pw.close()
    spark.stop()
   }

  
  def main(args:Array[String]) = {
      OParser.parse(arg_parser, args, Config()) match {
      case Some(config) => {
        phenotype_recommendation(config)
      }
      case _ => ()
    }
  }
  }
