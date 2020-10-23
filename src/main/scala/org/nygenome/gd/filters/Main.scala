package org.nygenome.gd.filters
import scala.language.reflectiveCalls
import org.rogach.scallop._

/*
    Filters Argument Specification

    TODO: better documentation of ..

    1) filter rules
    2) explanation of parameterization
*/

class Conf(arguments: Seq[String]) extends ScallopConf(arguments) {
  val filter_1 = new Subcommand("f1") {
    val tsv = opt[String](
      name="tsv",
      required = true,
      descr = "TSV (variants) from Exomiser output")

    val filter_config = opt[String](
      name = "filter_config",
      required = true,
      descr = "Filter configuration file (YAML) ")

    val input = opt[String](
      name = "input",
      required = true,
      descr = "Input rare VCF file from Exomiser output")

    val output = opt[String](
      name = "output",
      required = true,
      descr = "Output VCF after Filtering")
  }
  val filter_2 = new Subcommand("f2") {
    val tsv = opt[String](
      name="tsv",
      required = true,
      descr = "TSV (genes) from Exomiser output")

    val filter_config = opt[String](
      name = "filter_config",
      required = true,
      descr = "Filter configuration file (YAML) ")

    val exomiser_config = opt[String](
      name = "exomiser_config",
      required = true,
      descr = "Exomiser configuration file (YAML)")

    val sex = opt[String](
      name = "sex",
      required = true,
      descr = "Reported Sex")

    val input = opt[String](
      name = "input",
      required = true,
      descr = "Input variant prioritized VCF file")

    val output_csv = opt[String](
      name = "output_csv",
      required = true,
      descr = "Output [CSV]")

    val output_txt = opt[String](
      name = "output_txt",
      required = true,
      descr = "Output [TXT]")
  }
  val filter_3 = new Subcommand("f3") {
    val filter_config = opt[String](
      name = "filter_config",
      required = true,
      descr = "Filter configuration file (YAML) ")

    val tsv = opt[String](
      name="tsv",
      required = true,
      descr = "TSV (genes) from Exomiser 2 output")

    val parent_tsv = opt[String](
      name="parent_tsv",
      required = true,
      descr = "TSV (genes) from Exomiser 1 output")

    val input_vcf = opt[String](
      name = "input_vcf",
      required = true,
      descr = "Input rare-refined VCF file from Exomiser output")

    val output_vcf = opt[String](
      name = "output_vcf",
      required = true,
      descr = "Output Variants [VCF]")

    val output_csv = opt[String](
      name = "output_csv",
      required = true,
      descr = "Output Variants [CSV]")
  }

  addSubcommand(filter_1)
  addSubcommand(filter_2)
  addSubcommand(filter_3)

  verify()
}

object Main extends App {
  val conf = new Conf(args)
  conf.args.head match {
    case "f1" => {
      val c = conf.filter_1
      VariantFilter.main(
        c.filter_config(),
        c.tsv(),
        c.input(),
        c.output())
    }
    case "f2" => {
      val c = conf.filter_2
      PhenotypeFilter.main(
        c.filter_config(),
        c.exomiser_config(),
        c.tsv(),
        c.sex(),
        c.input(),
        c.output_csv(),
        c.output_txt()
      )}
    case "f3" => {
      val c = conf.filter_3
      ExtractionFilter.main(
        c.filter_config(),
        c.tsv(),
        c.parent_tsv(),
        c.input_vcf(),
        c.output_vcf(),
        c.output_csv()
      )
    }
  }
}
