name := "gd-filters"
version := "0.1.1-beta"
description := "Application to prioritize variants using phenotype information (HPO)"
maintainer:= "Kevin Shi <kshi@nygenome.org>"
packageSummary := "Genome Diver"
packageDescription := description.value
scalaVersion := "2.12.8"

scalacOptions ++= Seq("-unchecked", "-deprecation", "-feature", "-target:jvm-1.8")
javacOptions  ++= Seq("-source", "1.8", "-target", "1.8", "-Xlint")

// force JDK to 1.8, anything else has odd dependency issues
initialize := {
  val _ = initialize.value
  val required = "1.8"
  val current  = sys.props("java.specification.version")
  assert(current == required, s"Unsupported JDK: java.specification.version $current != $required")
}

// official scala repositories
resolvers += "Sonatype Releases" at "https://oss.sonatype.org/content/repositories/releases/"
resolvers += "Sonatype Snapshots" at "https://oss.sonatype.org/content/repositories/snapshots/"
resolvers += Resolver.sonatypeRepo("snapshots")

// custom libraries (external jars)
unmanagedBase := baseDirectory.value / "lib"

// sbt packaging to rpm, docker ...
enablePlugins(JavaServerAppPackaging)
enablePlugins(LinuxPlugin)

scalacOptions += "-Ypartial-unification"

libraryDependencies ++= Seq(
  "com.mchange"                  %% "danburkert-continuum"   % "0.3.99",
  "org.rogach"                   %% "scallop"                % "3.3.1",
  "org.typelevel"                %% "cats-core"              % "1.1.0",
  "io.circe"                     %% "circe-core"             % "0.10.0",
  "com.typesafe.scala-logging"   %% "scala-logging"          % "3.9.2",
  "io.circe"                     %% "circe-core"             % "0.10.0",
  "com.github.samtools"          %  "htsjdk"                 % "2.18.2",
  "io.circe"                     %% "circe-optics"           % "0.11.0",
  "io.circe"                     %% "circe-generic"          % "0.11.0",
  "io.circe"                     %% "circe-parser"           % "0.11.0",
  "io.circe"                     %% "circe-yaml"             % "0.9.0",
  "io.github.pityka"             %% "saddle-core-fork"       % "1.3.4-fork1",
  "org.monarchinitiative.phenol" % "phenol-core"             % "1.3.3",
  "org.monarchinitiative.phenol" % "phenol-io"               % "1.3.3"
)

cancelable in Global := true



