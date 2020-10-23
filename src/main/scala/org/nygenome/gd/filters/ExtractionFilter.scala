package org.nygenome.gd.filters

// Language, Logging
import scala.language.postfixOps
import scala.language.reflectiveCalls
import scala.collection.JavaConverters._
import com.typesafe.scalalogging.Logger
import cats.implicits._

// IO (Java), VCF Reader
import java.io.{File, FileInputStream, FileNotFoundException, InputStreamReader}
import htsjdk.tribble.AbstractFeatureReader
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.Genotype
import htsjdk.variant.variantcontext.writer.Options
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.vcf.VCFCodec
import htsjdk.variant.vcf.VCFHeader

// Config parser [YAML]
import io.circe.optics.JsonPath._
import io.circe.yaml.parser

// DataFrame
import org.saddle._
import org.saddle.io._
import org.saddle.io.CsvImplicits._

// Interval
import continuum.{Interval, IntervalSet}

object ExtractionFilter {
  val logger = Logger("[Variant Extraction]")

  // Extraction Filter
  // ---------------------------------------------------------------------------------
  // After Exomiser is ran again with updated inputs (VCF [rare], HPO [caregiver updated])
  //
  def main(_filter_config: String, _genes_tsv: String, _parent_genes_tsv:String, _input_vcf: String, _output_vcf: String,
           _output_csv: String): Unit = {

    logger.info("loading configurations")

    val YAML_FILTER_CONFIG = new File(_filter_config)
    val PARENT_GENES_TSV = new File(_parent_genes_tsv)
    val GENES_TSV = new File(_genes_tsv)
    val INPUT_VCF = new File(_input_vcf)
    val OUTPUT_CSV = new File(_output_csv)
    val OUTPUT_VCF = new File(_output_vcf)
    val FILTER_CONF = parser.parse(new InputStreamReader(new FileInputStream(YAML_FILTER_CONFIG))).right.get

    // error, sanity checking ...
    if (!GENES_TSV.exists()) {
      throw new FileNotFoundException("Genes TSV file missing")
    }
    if (!YAML_FILTER_CONFIG.exists()) {
      throw new FileNotFoundException("YAML filter file missing")
    }
    if (!INPUT_VCF.exists()) {
      throw new FileNotFoundException("Input VCF missing")
    }

    val GENES_TSV_PARSED = CsvFile(GENES_TSV.getAbsolutePath)
    //val GENES_EVIDENCE = List(0, 12, 13, 14, 15, 16, 17)
    val GENES_EVIDENCE = List(0, 12)

    // Common Constants / Keys
    val DISEASES = "DISEASES"
    val GENE = "GENE"
    val POS = "POS"
    val EXOMISER_VARIANT_SCORE = "EXOMISER_VARIANT_SCORE"
    val EXOMISER_GENE_PHENO_SCORE = "EXOMISER_GENE_PHENO_SCORE"
    val COMBINED_SCORE = "COMBINED_SCORE"
    val SCORE_MASK = Vec(true, true, true, true, true, true, true, true, true, false, false)

    val OMIM_ID = "OMIM_ID"
    val OMIM_REGEX = """(\d+)""".r
    val DISEASE_REGEX = """((?:ORPHA|OMIM):\d+)""".r
    val TOP_V_VARIANTS = root.filter_3.top_v.int.getOption(FILTER_CONF).get

    def DESC[T: Ordering] = implicitly[Ordering[T]].reverse

    logger.info("Loading a VCF")
    val vcf_reader_input = AbstractFeatureReader.getFeatureReader((INPUT_VCF).getAbsolutePath(), new VCFCodec(), false)
    val v = vcf_reader_input.iterator.asInstanceOf[java.util.Iterator[VariantContext]].asScala.foldLeft(((
      Vec("#CHROM"), Vec("POS"), Vec("REF"), Vec("ALT"), Vec("HGVS"), Vec("ZYGOSITY"),
      Vec("VARIANT_EFFECT"), Vec("OMIM_ID"), Vec("GENE_PHENO_SCORE"), Vec("GENE"), Vec("COMBINED_SCORE"),
    ),
      collection.mutable.Map[String, IntervalSet[Int]]().empty
    )) {
      (acc, vc) => {

        // exclude variants that are FILTERED
        // exclude variants that are marked "ExWarn"
        if (vc.isNotFiltered &&
          !vc.hasAttribute("ExWarn") && vc.hasAttribute("ExGeneSymbol") &&
          vc.hasAttribute("ExVarHgvs") && vc.hasAttribute("ExVarScore") &&
          vc.hasAttribute("ExGeneSPheno") && vc.hasAttribute("ExGeneSCombi")) {

          vc.getGenotypes().iterator.asInstanceOf[java.util.Iterator[Genotype]].asScala.foldLeft(acc) { (acc_gc, gc) => {
            val ref = vc.getReference.getBaseString
            val unique_alts = gc.getAlleles().asScala
              .filter(_.isNonReference)
              .map(_.getBaseString).distinct

            unique_alts.foldLeft(acc_gc) { (acc_alt, alt) => {
              var range = acc_alt._2 // mutable map, passed in the acc
              val accu = acc_alt._1 // Vector tuple above...
              val gene = vc.getAttribute("ExGeneSymbol").toString
              val pheno_score = vc.getAttribute("ExGeneSPheno").toString
              val ex_final_score = vc.getAttribute("ExGeneSCombi").toString
              val variant_idx = unique_alts.indexOf(alt)

              val variant_score = vc.getAttribute("ExVarScore").toString match {
                case arr if arr.startsWith("[") => arr.replace("[", "").replace("]", "").split(",")(variant_idx)
                case str => str
              }

              // rely on Exomiser Combined Score 
              val combined = ex_final_score.toFloat
              //val combined = (pheno_score.toFloat + variant_score.toFloat) / 2

              // essentially construct an set of variant diffs
              val variant_range = (vc.getStart to (vc.getStart + Math.max(ref.length - 1, alt.length - 1)))
              range(gene) = range.get(gene) match {
                case Some(s) => s + Interval.fromRange(variant_range)
                case None => IntervalSet(Interval.fromRange(variant_range))
              }

              val omim = (vc.hasAttribute("OMIM_MENDELIAN_MUT_ID") match {
                case true => "OMIM:" + OMIM_REGEX.findFirstIn(
                  vc.getAttribute("OMIM_MENDELIAN_MUT_ID").toString).mkString
                case false => ""
              })

              val hgvs = vc.getAttribute("ExVarHgvs").toString match {
                case arr if arr.startsWith("[") =>
                  arr.replace("[", "").replace("]", "").split(",")(variant_idx)
                case str => str
              }

              val var_eff = vc.getAttribute("ExVarEff").toString match {
                case arr if arr.startsWith("[") =>
                  arr.replace("[", "").replace("]", "").split(",")(variant_idx)
                case str => str
              }

              ((accu._1 concat Vec(vc.getContig),
                accu._2 concat Vec(vc.getStart.toString),
                accu._3 concat Vec(vc.getReference.getBaseString),
                accu._4 concat Vec(alt),
                accu._5 concat Vec(hgvs),
                accu._6 concat Vec(gc.getType().toString),
                accu._7 concat Vec(var_eff),
                accu._8 concat Vec(omim),
                accu._9 concat Vec(pheno_score.toString),
                accu._10 concat Vec(gene),
                accu._11 concat Vec(combined.toString)),
                range)
            }
            }
          }
          }
        } else {
          acc // filtered out 0
        }
      }
    }
    vcf_reader_input.close

    logger.info("Converting to a Dataframe")
    val frame = v._1

    // intervals represent a map of "overlapping variants"
    val intervals = v._2
    val v_df = Frame(
      frame._1, frame._2, frame._3, frame._4, frame._5,
      frame._6, frame._7, frame._8, frame._9, frame._10, frame._11)
      .withColIndex(0)
      .sortedRowsBy(r => r.get(COMBINED_SCORE).get.toFloat)(DESC)

    // top scoring genes (w/ minimum variant score cut off)
    val gene_top_combined = v_df
      .rmask(SCORE_MASK).squeeze.mapValues(_.toString)
      .withRowIndex(0).mapValues(_.toFloat)
      .groupBy.combine(_.toSeq.take(2).min)

    // mutable map containing variants as gene -> IntervalSet; used to coalesce variants
    // that should be a singular call
    var included_intervals = collection.mutable.Map[String, IntervalSet[Int]]().empty

    def update_included(gene: String, reference_interval: IntervalSet[Int], pos: Int): Unit = {
      synchronized {
        reference_interval.unioning(Interval.point(pos)).headOption match {
          case Some(i) => included_intervals.get(gene) match {
            case Some(inc_set) => included_intervals(gene) = inc_set + i
            case None => included_intervals(gene) = IntervalSet(i)
          }
          case None => // do nothing, optimistic update.
        }
      }
    }

    val top_variants_vcf = v_df
      .rfilter(r => {
        val gene         = r.get(GENE).getOrElse("-").toString
        val pos          = r.get(POS).getOrElse(0).toString.toInt
        val score        = r.get(COMBINED_SCORE).getOrElse("0").toString.toFloat
        val min_score    = gene_top_combined.row(gene)

        // check if minimum score exist and check if intervals exists for gene
        (min_score.isEmpty, intervals.get(gene)) match {
          case (false, Some(variant_interval)) => {
            (score >= min_score.raw(0,0).toFloat) match {
              case true => {
                included_intervals.get(gene) match {
                  case Some(included_interval_set) => {
                    included_interval_set.contains(Interval.point(pos)) match {

                      // if variant is not in included then potentially add its interval to included interval
                      case false => {
                        if (included_interval_set.size < 2) {
                          update_included(gene, variant_interval, pos)
                          true
                        } else {
                          // max number of variant recommendations  reached
                          false
                        }
                      }

                      // if variant is in included set then filter out
                      // already established the canonical variant representing
                      // the grouping
                      case true => false
                    }
                  }
                  case None => {
                    // if included doesn't exists for gene add gene and
                    // the first interval
                    update_included(gene, variant_interval, pos)
                    true
                  }
                }
              }
              // minimum score not met
              case _ => false
            }
          }
          // not enough information (min_score missing)
          case (_, _)  => false
        }
      }).rowSlice(0, TOP_V_VARIANTS)

    // Find relevant diseases (according to Exomiser)
    // load all evidence in genes.tsv
    logger.info("Extracting relevant diseases")
    val relevant_diseases =
      CsvParser.parse(GENES_EVIDENCE, CsvParams('\t'))(GENES_TSV_PARSED)
        .withRowIndex(0).withColIndex(0)
        .mapValues(v => {
          // Extract out Normalized Disease Terms
          DISEASE_REGEX.findAllIn(v.toString).mkString(",")
        })
        .rfilterIx {
          // Filter for implicated genes (from variants file)
          top_variants_vcf
            .rmask(SCORE_MASK).squeeze.mapValues(_.toString)
            .withRowIndex(0).rowIx.contains(_)
        }
        .rreduce(r => {
          // Flatten the columns
          r.toVec.toSeq.map(_.trim).toSet.filter(!_.isEmpty).mkString(", ")
        })

    val completion = top_variants_vcf
      .rtransform(r => {
        val new_idx = Index(DISEASES)
        val gene = r.get(GENE).get.toString
        val disease = relevant_diseases.row(gene)

        r.concat(disease.isEmpty match {
          case false => Series(Vec(disease.raw(0, 0)), new_idx)
          case true => Series(Vec(""), new_idx)
        })
      })

    // Ignore #CHROM, POS, REF, ALT for CSV purposes
    // - need those values only for VCF filtering


    // TODO.. output do score diff here:
    def csv_diff_df(csv_file:File):Frame[String, String, Double] = {
      CsvParser.parse(List(0, 4), CsvParams('\t'))(CsvFile(csv_file.getAbsolutePath))
        .withRowIndex(0).withColIndex(0)
        .mapValues(_.toDouble)
    }

    val diff_merged  = csv_diff_df(PARENT_GENES_TSV)
      .joinPreserveColIx(csv_diff_df(GENES_TSV), how=index.InnerJoin)
      .toRowSeq.foldLeft(Frame.empty[String, String, Double]){
        case (frames, element) =>
          val (key, series) = element
          val next = Frame(key ->
            Series(Vec(series.at(1) - series.at(0)),
              Index("delta")))
          frames.rconcat(next)
    }.T

    val output_csv = completion
      .rmask(Vec(true, true, true, true, false, false,
        false, true, false, false, false, false))
      .squeeze
      .rtransform(r => {
        val new_idx = Index("SCORE_CHANGE")
        val gene = r.get(GENE).get.toString
        val score_delta = diff_merged.row(gene)

        r.concat(score_delta.isEmpty match {
          case false => Series(Vec(score_delta.raw(0, 0).toString), new_idx)
          case true => Series(Vec(""), new_idx)
        })
      })

    // ~ flush (CSV)
    logger.info("Output variants as CSV")
    output_csv.writeCsvFile(OUTPUT_CSV.getAbsolutePath, withColIx = true, withRowIx = true)

    logger.info("Output variants as VCF")
    val vcf_reader = AbstractFeatureReader.getFeatureReader(
      INPUT_VCF.getAbsolutePath(),
      new VCFCodec(), false)

    val vcf_writer = new VariantContextWriterBuilder()
      .setOutputFile(OUTPUT_VCF)
      .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
      .unsetOption(Options.INDEX_ON_THE_FLY).build();

    vcf_writer.writeHeader(vcf_reader.getHeader.asInstanceOf[VCFHeader])
    for (vc: VariantContext <- vcf_reader.iterator.asInstanceOf[java.util.Iterator[VariantContext]].asScala) {
      for (alt <- vc.getAlternateAlleles.asScala.map(_.getBaseString)) {
        vc.isNotFiltered && !completion.rfilter(r => {
          r.get("#CHROM").getOrElse("") == vc.getContig &&
            r.get("POS").getOrElse("") == vc.getStart.toString &&
            r.get("REF").getOrElse("") == vc.getReference.getBaseString &&
            r.get("ALT").getOrElse("") == alt

        }).isEmpty match {
          case true => vcf_writer.add(vc)
          case _ => {
            // filtered out
          }
        }
      }
    }

    // ~ flush (VCF)
    vcf_writer.close()
  }
}
