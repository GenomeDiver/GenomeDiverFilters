package org.nygenome.gd.filters

// Language, Logging
import java.io.PrintWriter

import scala.language.postfixOps
import scala.language.reflectiveCalls
import scala.collection.JavaConverters._
import scala.collection.{Set, SortedSet}
import scala.collection.mutable.{Set => MutableSet}
import com.typesafe.scalalogging.Logger
import cats.implicits._

// IO (Java), VCF Reader
import java.io.{File, FileInputStream, FileNotFoundException, InputStreamReader}
import htsjdk.tribble.AbstractFeatureReader
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.{VCFCodec}

// Config parser [YAML]
// import io.circe._
import io.circe.optics.JsonPath._
import io.circe.yaml.parser

// Ontology [HPO]
import org.monarchinitiative.phenol.io.OntologyLoader
import org.monarchinitiative.phenol.io.obo.hpo.HpoDiseaseAnnotationParser
import org.monarchinitiative.phenol.ontology.algo.OntologyAlgorithm.{getDescendents => OntologyGetDescendents}
import org.monarchinitiative.phenol.ontology.data.TermId
import org.monarchinitiative.phenol.formats.hpo.HpoAnnotation

// DataFrame
import org.saddle._
import org.saddle.io._
import org.saddle.io.CsvImplicits._
import org.saddle.scalar.{Scalar, NA}

object PhenotypeFilter {

  val logger = Logger("[Phenotype Prioritization]")

  // Phenotype Prioritization Filter
  // ---------------------------------------------------------------------------------
  //

  // Configurations
  // ---------------------------------------------------------------------------------
  // Files:
  //
  //    Genes TSV            Annotation of implicated genes (directly from Exomiser),
  //                         contains gene phenotype score.
  //    Input VCF:           Filtered VCF (1) provides the relevant variants
  //    Filer YAML:          File for constants and paths to HPO ontology/annotations
  //    Exomiser YAML:       File (used to isolate the original phenotype entered)
  //    Output CSV:          List of Phenotypes as determined by Categorization 1 + 2
  //
  // Parameters (from filters.yml)
  //
  //    REPORTED SEX         FEMALE | MALE | UNKNOWN
  //    TOP_K:               Top K unique variants scores
  //    MAX_PHENO:           The maximum number of phenotypes to show
  //    INPUT_HPO:           List of HPO terms from the Exomiser Input
  //    VSCORE_MIN_CAT1:     Variant score cutoff for Categorization 1
  //    VSCORE_MIN_CAT2:     Variant score cutoff for Categorization 2
  //    HPO_ONTOLOGY_PATH    Path to the reference HPO ontology (.obo)
  //    HPO_ANNOTATION_PATH: Path to the reference HPO (gene - phenotype) annotation
  // ---------------------------------------------------------------------------------
  def main(filter_config:String, exomiser_config:String, tsv:String, sex:String, input:String, output_csv:String, output_txt:String): Unit = {

    logger.debug("loading configurations")

    val YAML_FILTER_CONFIG = new File(filter_config)
    val YAML_EXOMISER_CONFIG = new File(exomiser_config)
    val INPUT_VCF = new File(input)
    val OUTPUT_CSV = new File(output_csv)
    val OUTPUT_TXT = new File(output_txt)

    // error/sanity checking ...
    if (!(new File(tsv)).exists()) {
      throw new FileNotFoundException("TSV file missing")
    }
    if (!YAML_FILTER_CONFIG.exists()) {
      throw new FileNotFoundException("YAML filter file missing")
    }
    if (!YAML_EXOMISER_CONFIG.exists()) {
      throw new FileNotFoundException("YAML exomiser file missing")
    }
    if (!INPUT_VCF.exists()) {
      throw new FileNotFoundException("Input VCF missing")
    }

    val TSV_GENES = CsvFile(tsv)
    val FILTER_CONF = parser.parse(new InputStreamReader(new FileInputStream(YAML_FILTER_CONFIG))).right.get
    val EXOMISER_CONF = parser.parse(new InputStreamReader(new FileInputStream(YAML_EXOMISER_CONFIG))).right.get

    val HPO_ONTOLOGY_PATH = root.filter_2.hpo_ontology.string.getOption(FILTER_CONF).get
    val HPO_ANNOTATION_PATH = root.filter_2.hpo_annotation.string.getOption(FILTER_CONF).get
    val HPO_ANNOTATION_DISEASE_PATH = root.filter_2.hpo_annotation_disease.string.getOption(FILTER_CONF).get
    val TOP_K: Int = root.filter_2.top_k.int.getOption(FILTER_CONF).get
    val VSCORE_MIN_CAT1: Double = root.filter_2.variant_score_min_cat_1.double.getOption(FILTER_CONF).get
    val VSCORE_MIN_CAT2: Double = root.filter_2.variant_score_min_cat_2.double.getOption(FILTER_CONF).get
    val MAX_PHENO: Int = root.filter_2.max_pheno.int.getOption(FILTER_CONF).get
    val INPUT_HPO: List[String] = root.analysis.hpoIds.each.string.getAll(EXOMISER_CONF)
    val PHENOTYPE_ABNORMALITY = "HP:0000118"
    val REPORTED_SEX = sex

    // Sort Ordering
    def DESC[T: Ordering] = implicitly[Ordering[T]].reverse

    logger.debug("Loading HPO ontology")
    val ontology = OntologyLoader.loadOntology(new File(HPO_ONTOLOGY_PATH))

    logger.debug("Loading HPO annotations (disease -> sex/phenotype) + normalizing input")
    val phenotype_to_sex = CsvParser.parse(List(3, 8), CsvParams(separChar = '\t'))(CsvFile(HPO_ANNOTATION_DISEASE_PATH))
      .withColIndex(0).withRowIndex(0)
      .rfilter { case r => !r.at(0).get.trim.isEmpty }
      .mapValues(_.toUpperCase)

    logger.debug("Loading HPO annotations (disease -> phenotype)")
    val disease = new HpoDiseaseAnnotationParser(HPO_ANNOTATION_DISEASE_PATH, ontology)
    val disease_info = disease.parse()
    val hpo_to_disease_term = disease.getTermToDiseaseMap()
    val total_annotated_diseases = CsvParser.parse(List(0, 1), CsvParams(separChar = '\t'))(CsvFile(HPO_ANNOTATION_DISEASE_PATH))
      .withColIndex(0).withRowIndex(0).rowIx.uniques.length.toDouble

    logger.debug("Loading HPO annotations (gene -> phenotype)")
    val gene_pheno_ref = CsvParser.parse(List(1, 3), CsvParams(separChar = '\t'))(CsvFile(HPO_ANNOTATION_PATH))
      .setColIndex(Index("GENE_SYMBOL", "HPO_ID"))
      .withRowIndex(0)

    logger.debug("Loading in genes, disease recommendations from exomiser (genes.tsv)")
    val df_g_ref = CsvParser.parse(List(0, 2, 12, 13, 14, 15, 16, 17), CsvParams('\t'))(TSV_GENES)

    logger.debug("Extracting diseases (ORPHA, OMIM) from evidence file - HUMAN_PHENO evidence")
    val total_evidence = df_g_ref.where(Vec(true, false, true, false, false, false, false, false))
      .withRowIndex(0).withColIndex(0)
      .mapValues(v => {
        """((?:ORPHA|OMIM):\d+)""".r.findAllIn(v.toString).mkString(",")
      })

    // ---------------------------------------------------------------------------------
    // Compute "distinct phenotypes" across genes
    //
    // Example:
    //    geneA: [HP:001, HP:002, ...]
    //    geneA: [(HP:001, 1), (HP:002, 20), ...]
    //    geneA: [(HP:002, 20), (HP:001, 1), ...]
    //    geneA: HP:002
    //
    // ---------------------------------------------------------------------------------
    def distinct_phenotypes(gene_pheno: Frame[String, String, Any], blacklist: Set[String] = Set.empty[String]):
    Frame[String, String, Scalar[String]] = {

      // Set of all descendants for a list of HPO IDs
      def descendants(hpo: Iterable[String]): Set[String] = {
        hpo.foldLeft(Set.empty[String]) { (acc, hpo_id) => {
          acc ++ OntologyGetDescendents(ontology, TermId.of(hpo_id)).asScala.toSeq.map(_.toString).toSet
        }
        }
      }

      def ascendants(hpo: Iterable[String]): Set[String] = {
        hpo.foldLeft(Set.empty[String]) { (acc, hpo_id) => {
          acc ++ ontology.getAncestorTermIds(TermId.of(hpo_id)).asScala.toSeq.map(_.toString).toSet
        }
        }
      }

      def toDefaultPrecision(value: BigDecimal): Double = {
        value.setScale(4, BigDecimal.RoundingMode.HALF_UP).toDouble
      }

      // some convenience classes
      case class HpoScore(hpo: String, score: Double)
      case class GenePheno(gene: String, pheno: String)

      // find only diseases affiliated with interested genes
      val valid_diseases = total_evidence.rfilterIx {
        gene_pheno.rowIx.contains(_)
      }
        .toMat.rows.foldLeft(Set.empty[String]) {
        (acc, v) => {
          if (!v.isEmpty) {
            acc ++ v.toSeq.map(_.trim).toSet.filter(!_.isEmpty)
          } else acc
        }
      }

      logger.info(s"valid diseases: ${valid_diseases.toString}")
      logger.info("\nHPO_ID\tIC\tDIS\tfDIS\tSUM_FREQ|fDIS\t(IC*(SUM_FREQ|fDIS) [GD PHENO SCORE])")

      // [GENE -> HPO_ID, GENE_PHENO_SCORE, EXOMISER_VARIANT_SCORE]
      // =====================================================================
      // - filter out phenotypes in the input set
      // - filter out and phenotypes in the ascendant chain of the input set
      // - filter out phenotypes not considered "phenotypic abnormality"
      // - filter out phenotypes when sex reported and annotated
      // - calculate the GD score to all recommended phenotypes
      // - attribute the GD score to all recommended phenotypes
      val scored_phenotypes = gene_pheno.map {
        case (r, c, v) =>
          c match {
            case "HPO_ID" =>
              (r, c, v.toString.split(",").map(_.trim)
                .filter(hpo => !(INPUT_HPO.toSet.contains(hpo) || ascendants(INPUT_HPO).contains(hpo.toString)))
                .filter(hpo => ascendants(List(hpo)).contains(PHENOTYPE_ABNORMALITY))
                //.filter(hpo => {
                // if EITHER reported or annotations do not exist, THEN pass
                // if BOTH reported AND annotations exist AND disagree then reject phenotype
                // TODO: intersex situation? discordant sex as a flag?
                //                  (REPORTED_SEX, phenotype_to_sex.first(hpo).get("Sex")) match {
                //                    case ("UNKNOWN", _)               => true
                //                    case (_, NA)                      => true
                //                    case (x:String, y:Scalar[String]) => (x == y.get)
                //                  }
                //  true
                //})
                .map(hpo => {
                  val hpo_term = TermId.of(hpo)
                  val diseases = hpo_to_disease_term.get(hpo_term)
                  val filtered_diseases: List[TermId] = diseases.iterator.asScala.toList.filter(d => {
                    valid_diseases.contains(d.toString)
                  })

                  // phenotype frequency in disease.. usually phenotype exists in one disease.
                  // ex. [.5] or [.3, .8]
                  val freq = filtered_diseases.foldLeft(Set.empty[Double]) { (acc, x) => {
                    disease_info.get(x).getAnnotation(hpo_term) match {
                      case null => acc
                      case x: HpoAnnotation => acc + x.getFrequency()
                    }
                  }}

                  // deprecated
                  val sum_freq = toDefaultPrecision(BigDecimal(freq.sum))
                  val ic_num = if (diseases.size == 0) total_annotated_diseases else diseases.size.toDouble
                  val ic = toDefaultPrecision(BigDecimal(-Math.log(ic_num / total_annotated_diseases)))

                  // phenotype has some affinity to some of the interested diseases.
                  // take the average of these prevalence.
                  val avg_freq = freq.size match {
                    case 0 => toDefaultPrecision(BigDecimal(0.0))
                    case _ => toDefaultPrecision(BigDecimal(freq.sum/freq.size))
                  }

                  if (avg_freq > 0) {
                    logger.info(s"${hpo}\t${ic}\t${diseases.size}\t${filtered_diseases.size}\t${sum_freq}\t${avg_freq}")
                  }

                  // attribute [GD PHENO SCORE] to HPO phenotype
                  HpoScore(hpo, avg_freq)
                }).sortBy {
                _.score
              }(DESC)
              )
            case _ => (r, c, v)
          }
      }

      // Draft Picks of Phenotypes
      // =====================================================================
      // 1st draft - order genes by highest GD pheno score   => submit high scoring phenotypes to draft
      // 2nd draft - order genes by variant score            => submit high scoring phenotypes to draft
      var gene_pheno_set = MutableSet.empty[GenePheno]

      // Possibility of an empty variant results due to filtering
      if (!scored_phenotypes.isEmpty) {

        // 1st Draft Pick
        // ---------------------------------------------------------------------
        // 1) order by the highest (GD PHENO SCORE)
        val first_draft = scored_phenotypes
          .sortedRowsBy(r => {
            val scores_at_gene = r.at(0).get.asInstanceOf[Array[HpoScore]]
            scores_at_gene.headOption.getOrElse(0).asInstanceOf[HpoScore].score
          })(DESC).colAt(loc = 0).toSeq

        // 2) submit highest scoring phenotype into draft
        //   * make sure phenotypes cannot be selected twice *
        first_draft.foreach { case (gene, pheno_list) => {
          pheno_list.asInstanceOf[Array[HpoScore]]
            .filter(hs => !blacklist.contains(hs.hpo))
            .headOption match {
              case Some(a) => {
                if (!gene_pheno_set.map(_.pheno).contains(a.hpo)) {
                  synchronized {
                    gene_pheno_set.add(GenePheno(gene, a.hpo))
                  }
                }
              }
              case None =>
              // phenotypes exhausted for gene
              // hopefully picked up by nth-draft
            }
        }}

        // (N+1) Draft Pick
        // ---------------------------------------------------------------------
        // 1) Order by variant score
        // 2) Filter out genes from the first pick
        // 3) Find the highest GD Pheno Score
        val nth_draft = scored_phenotypes
          .sortedRowsBy(r => {
            r.at(loc = 2).get.asInstanceOf[Float]
          })(DESC)
          // .rfilterIx {!gene_pheno_set.map(_.gene).contains(_)}
          .colAt(0).toSeq

        // run draft x2
        1 to 4 foreach { _ =>
          nth_draft.foreach { case (gene, pheno_list) => {
            pheno_list.asInstanceOf[Array[HpoScore]]
              .filter(hs => !gene_pheno_set.map(_.pheno).contains(hs.hpo))
              .filter(hs => !blacklist.contains(hs.hpo))
              .headOption match {
              case Some(a) => {
                synchronized {
                  gene_pheno_set.add(GenePheno(gene, a.hpo))
                }
              }
              case None =>
              // all phenotypes exhausted for gene
              // ~ gene has no phenotype contribution to current categorization,
            }
          }}
        }

      }
      // Construct a Dataframe Result Object
      // =====================================================================
      // - transform the Set(GenePheno) into a DataFrame (gene -> recommendation [HPO])
      // - set "GENE" as index and make sure values are Scalar[String]
      //
      {
        val seq = gene_pheno_set.toArray
        Frame(
          "GENE" -> Vec(seq.map(v => v.gene): _*),
          "RECOMMENDATION" -> Vec(seq.map(v => v.pheno): _*))
          .withRowIndex(col = 0)
          .mapValues(v => Scalar(v))
      }
    }

    // ---------------------------------------------------------------------------------
    // CATEGORIZATION TIER (1)
    // ---------------------------------------------------------------------------------
    // Prioritize Genes that have a high connection to the disease
    //
    // Pick genes from Exomiser with any of the top (k)
    // gene_phenotype_score values
    //
    // [Filter]       Gene has variants >=1 variant in VCF-rare
    // [Filter]       If any of the variant of a gene has a variant score of 0.9
    // [Annotate]     For each shortlisted gene, find all phenotypes associated with that gene
    // [Filter]       Isolate "distinctive phenotypes"
    //
    // ---------------------------------------------------------------------------------

    /*  Load in all candidate genes (by Exomiser) for inclusion into short list using the Genes TSV file fromq Exomiser */
    logger.debug("loading candidate genes (from Exomiser)")
    val df_g = df_g_ref
      .where(Vec(true, true, false, false, false, false, false, false))
      .withRowIndex(col = 0).withColIndex(0)
      .mapValues(CsvParser.parseFloat)
      .sortedRowsBy(r => r.at(loc = 0))(DESC)

    /* Establish a ranking of TOP_K scores using [EXOMISER_GENE_PHENO_SCORE]
        Slice out the [TOP K (SCORED)] genes */
    val scores = df_g.colAt(loc = 0).toVec
    val rank = scores.foldLeftWhile(SortedSet.empty[Float])(_ + _)((s: SortedSet[Float], _) => s.size < TOP_K)
    val df_gs = df_g.rowSlice(from = 0, scores.findOne(v => v < rank.head))
    val sel_genes = df_gs.rowIx

    /* Load in variants of only the [TOP K (SCORED)] genes from [VCF RARE]
      1) MUST have PASS in the filter as [FILTER 1]
      2) MUST include relevant gene name (ExGeneSymbol)
      3) MUST include relevant variant score (ExVarScore)
      4) MUST not contain "ExWarn" tag

     i.e [GENE, VARIANT_SCORE]
     */
    logger.debug("loading (subset of) VCF file into a DataFrame")
    val vcf_reader = AbstractFeatureReader.getFeatureReader((INPUT_VCF).getAbsolutePath(), new VCFCodec(), false)
    val v = vcf_reader.iterator.asInstanceOf[java.util.Iterator[VariantContext]].asScala.foldLeft((Vec("EXOMISER_GENE"), Vec("EXOMISER_VARIANT_SCORE"))) {
      (acc, vc) => {

        // exclude variants that are FILTERED
        // exclude variants that are marked "ExWarn"
        if (vc.isNotFiltered && !vc.hasAttribute("ExWarn") &&
          vc.hasAttribute("ExGeneSymbol") && vc.hasAttribute("ExVarScore") &&
          sel_genes.contains(vc.getAttribute("ExGeneSymbol").toString)) {

          // include variants passing criteria...
          // exomiser sometimes gives TWO variants scores separated by comma, taking MAX score
          (acc._1 concat Vec(vc.getAttribute("ExGeneSymbol").toString),
            acc._2 concat Vec(vc.getAttribute("ExVarScore")
              .toString.replace("[", "").replace("]", "")
              .split(",").map(_.trim.toFloat).max.toString))
        } else {
          acc
        }
      }
    }
    vcf_reader.close

    // Group Variants by Genes and Scores aggregate on MAX([EXOME_VARIANT_SCORE]) of the [TOP K SCORED] genes
    // ex. [GENE (index), VARIANT_SCORE]
    val df_vg = Frame(v._1, v._2)
      .withColIndex(row = 0)
      .withRowIndex(col = 0).mapValues(CsvParser.parseFloat)
      .groupBy.combine(_.toSeq.reduceLeft(_ max _))

    // Filter out Genes where the variant score is beneath a certain MIN
    // ex. [GENE (index), PHENO_SCORE, VARIANT_SCORE]
    val df_scored_genes = df_gs.joinPreserveColIx(df_vg)
      .rfilter { case r => r.at(1) > VSCORE_MIN_CAT1 }

    // Get all the all the implicated genes and respective HPO phenotype annotations
    // ex. [GENE (index), HPO ID, PHENO_SCORE, VARIANT_SCORE]
    val df_gene_pheno = gene_pheno_ref
      .rfilterIx(df_scored_genes.rowIx.contains(_))
      .groupBy.combine(_.toSeq.mkString(","))

    logger.debug("calculating distinct phenotypes... [categorization 1]")

    // ex. [GENE (index), HPO ID, PHENO_SCORE, VARIANT_SCORE]
    val cat_1_all = distinct_phenotypes(df_gene_pheno.joinAnyPreserveColIx(df_scored_genes))

    // hide genes that have no phenotype recommendation
    val cat_1 = cat_1_all.joinAnyPreserveColIx(df_scored_genes, how = index.RightJoin)
      .rfilter(_.at(loc = 0).toString.contains(s"HP:"))
      .head(MAX_PHENO)

    // ---------------------------------------------------------------------------------
    // CATEGORIZATION TIER (2)
    // ---------------------------------------------------------------------------------
    // Prioritize Exomiser Variant Score (deleterious events in genes)
    //
    // 1) Categorization (2) only shows if Categorization (1) does not meet MAX_PHENO requirements
    // 2) Variant must pass a new threshold minimum (VSCORE_MIN_CAT2)
    // 3) Gene must not be in categorization tier (1)
    // 4) Find all associated phenotypes to gene subtracted by phenotypes in categorization tier (1)
    //
    //  * not all genes - have HPO annotation (interesting) hence the rfilter(!_.hasNA)
    // ---------------------------------------------------------------------------------
    // Determine if Categorization (2) is necessary
    val CAT2_MAX_PHENO: Int = MAX_PHENO - cat_1.numRows

    val result: Frame[String, String, Scalar[String]] = {
      if (CAT2_MAX_PHENO > 0) { // rule: 1

        logger.debug("calculating distinct phenotypes... [categorization 2]")

        // build categorization (2)
        val df_cat_2 = df_vg
          .rfilter {_.at(loc = 0) > VSCORE_MIN_CAT2}  // rule: 2
          .rfilterIx {!cat_1_all.rowIx.contains(_)}   // rule: 3
          .joinAnyPreserveColIx(gene_pheno_ref)       // acquire phenotypes (gene -> pheno)
          .where(Vec(false, true))                    // mask out scores
          .rfilter(!_.hasNA)                          // certain genes have no phenotypes? TODO: investigate
          .groupBy.combine(_.toSeq.mkString(","))     // merge phenotypes
          .mapValues(v => Scalar(v))

        // generate a blacklist of phenotypes from categorization (4)
        val cat_1_blacklist = cat_1_all.colAt(loc = 0).toSeq.map(_._2.toString).toSet
        val cat_2_scores = df_gs.joinPreserveColIx(df_vg).rfilter { case r => r.at(1) > VSCORE_MIN_CAT2 }
        val cat_2_all = distinct_phenotypes(df_cat_2.joinAnyPreserveColIx(cat_2_scores, how = index.LeftJoin), cat_1_blacklist)

        // hide genes that have no phenotype recommendation (taken by other genes)
        val cat_2 = cat_2_all.joinAnyPreserveColIx(cat_2_scores, how = index.RightJoin)
          .rfilter(_.at(0).toString.contains(s"HP:"))
          .head(CAT2_MAX_PHENO)

        // return categorization (1) and categorization (2) combined
        // * change type of scores (float) to string
        cat_1.concat(cat_2).mapValues(v => Scalar(v.toString))

      } else {

        // return sorted categorization (1) as quota is exceeded
        cat_1.mapValues(v => Scalar(v.toString))
      }
    }

    // write the resulting DataFrame to the output path
    result.writeCsvFile(OUTPUT_CSV.toPath.toString)

    {// output recommended HPOs to a file
      val pw = new PrintWriter(OUTPUT_TXT)
      pw.write(result.colAt(loc = 0).toSeq.foldLeft("") {
        case (acc, (_, rec: Scalar[String])) => acc + s"${rec.toString}\n"
      }.stripLineEnd)
      pw.close()
    }
  }
}

