package org.nygenome.gd.filters

import java.io.{File, FileInputStream, FileNotFoundException, InputStreamReader}

import scala.collection.JavaConverters._
import scala.language.postfixOps
import cats.implicits._
import com.typesafe.scalalogging.Logger
import htsjdk.tribble.AbstractFeatureReader
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.writer.{Options, VariantContextWriterBuilder}
import htsjdk.variant.vcf.{VCFCodec, VCFHeader}
import io.circe.optics.JsonPath.root
import io.circe.yaml.parser
import org.saddle.io.{CsvFile, CsvParams, CsvParser}

object VariantFilter {

  val logger = Logger("[Variant Prioritization]")

  // Variant Prioritization Filter
  // ---------------------------------------------------------------------------------
  //

  // Configurations
  // ---------------------------------------------------------------------------------
  // Files:
  //
  // Annotations TSV (used for REMM Score)
  // YAML config file contains externalized variables
  // Input VCF should be the Clinvar annotated Exomiser output (*.annotated.vcf)
  // Output VCF (filtered variants)

  def main(filter_config:String, tsv:String, input:String, output:String): Unit = {
    logger.debug("loading configurations")

    val TSV                             = new File(tsv)
    val YAML_CONFIG                     = new File(filter_config)
    val INPUT_VCF                       = new File(input)
    val OUTPUT_VCF                      = new File(output)
    val VARIANTS_TSV                    = CsvFile(tsv)
    val COLUMNS_TSV                     = List(0, 1, 2, 3, 5, 8, 15, 17)
    val FILTER_CONF                     = parser.parse(new InputStreamReader(new FileInputStream(YAML_CONFIG))).right.get
    val REMM_CUTOFF:Double              = root.filter_1.remm_min.double.getOption(FILTER_CONF).get
    val MAX_AF_CUTOFF:Double            = root.filter_1.af_max.double.getOption(FILTER_CONF).get
    val FUNCTIONAL_CLASSES:List[String] = root.filter_1.functional_classes.each.string.getAll(FILTER_CONF)

    if (!TSV.exists())                  {throw new FileNotFoundException("TSV file missing")}
    if (!YAML_CONFIG.exists())          {throw new FileNotFoundException("YAML file missing")}
    if (!INPUT_VCF.exists())            {throw new FileNotFoundException("Input VCF missing")}

    // ---------------------------------------------------------------------------------
    // Filter for
    //  - variants in the functional classes (from config)
    //  - [OR] variants that pass the REMM score threshold
    //
    // ALLELE FREQ AND (REMM > 0.5 OR (IN FUNC))
    // ---------------------------------------------------------------------------------
    logger.debug(s"extracting variants (REMM > ${REMM_CUTOFF}, MAX_FREQUENCY < ${MAX_AF_CUTOFF})")
    val passed_variants = CsvParser.parse(COLUMNS_TSV, CsvParams('\t'))(VARIANTS_TSV)
      .withColIndex(0)
      .rfilter{r =>

        // PASS FILTER AND [<MAX_AF AND [>REMM OR FUNCTION CLASS]]
        (!r.at(4).isNA && (r.at(4).get match {
          case "PASS" => true
          case _      => false
        })) && (
            (!r.at(7).isNA && (r.at(7).get match {
              case "."    => false
              case x      => x.toDouble <= MAX_AF_CUTOFF
            }))
            && (
            (!r.at(6).isNA && (r.at(6).get match {
              case "."    => false
              case x      => x.toDouble > REMM_CUTOFF
            })) ||
            (!r.at(5).isNA && (r.at(5).get match {
              case x => FUNCTIONAL_CLASSES contains (x toUpperCase)
            }))
          )
        )
      }
      .toMat.rows().foldLeft(
      Set.empty[Tuple4[String, String, String,String]]){
      (acc, v) => acc + Tuple4(v.raw(0), v.raw(1), v.raw(2), v.raw(3))
    }

    // ---------------------------------------------------------------------------------
    // Write to output VCF file
    // Before inclusion of variants to examine Clinvar assertion for a second opinion
    //  -  BENIGN                => discard
    //  -  BENIGN + CONFLICTING  => keep
    // ---------------------------------------------------------------------------------
    logger.debug("Writing VCF output")
    val vcf_reader = AbstractFeatureReader.getFeatureReader(
      INPUT_VCF.getAbsolutePath(), new VCFCodec(), false)

    val vcf_writer = new VariantContextWriterBuilder()
      .setOutputFile(OUTPUT_VCF)
      .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
      .unsetOption(Options.INDEX_ON_THE_FLY).build();

    vcf_writer.writeHeader(vcf_reader.getHeader.asInstanceOf[VCFHeader])
    for (vc:VariantContext <- vcf_reader.iterator.asInstanceOf[java.util.Iterator[VariantContext]].asScala) {
      Option(vc.getContig) match {
        case Some(chr) => {
          for (alt <- vc.getAlternateAlleles.asScala.map(_.getBaseString)) {
            val loc = Tuple4(chr, vc.getStart.toString, vc.getReference.getBaseString, alt)
            val clinvar_pathogenic = vc.getAttributeAsString("clinvar_pathogenic", "").toLowerCase
            val clinvar_review     = vc.getAttributeAsString("clinvar_review", "").toLowerCase

            // variant is FILTER:PASS
            if (vc.isNotFiltered) {

              (clinvar_pathogenic contains "pathogenic") match {

                // clinvar pathogenic OR likely pathogenic => add variant
                // - do not be concerned with review_status
                case true  => vcf_writer.add(vc)
                case false => {

                  // 1st filter, AF, REMM, FUNC classes (passed_variants
                  // 2nd filter, check if BENIGN
                  if (passed_variants.contains(loc)) {

                    val clinvar_combined_assertion = (
                      (clinvar_pathogenic equals "benign"),
                      (clinvar_review contains "conflicting")
                    )

                    clinvar_combined_assertion match {

                      // do nothing, clinvar confidently BENIGN             => do nothing
                      case (true, false)  =>

                      // clinvar BENIGN or LIKELY BENIGN, but CONFLICTING   => add variant
                      case (true, true)   => vcf_writer.add(vc)

                      // no clinvar BENIGN assertion                        => add variant
                      case (false, _)     => vcf_writer.add(vc)
                      case _ =>
                    }
                  }
                }
              }
            }
          }
        }
        case None => {/* do nothing - chromosome => None */}
      }
    }
    vcf_writer.close()
  }
}