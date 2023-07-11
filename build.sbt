import Dependencies._

ThisBuild / scalaVersion     := "2.13.10"
ThisBuild / version          := "0.1.0-SNAPSHOT"
ThisBuild / organization     := "org.nygenome"
ThisBuild / organizationName := "New York Genome Center"

lazy val root = (project in file("."))
  .settings(
    name := "recommendation",
    libraryDependencies += munit % Test
  )

unmanagedBase := baseDirectory.value / "libs"

libraryDependencies ++= Seq(
"dev.zio" %% "zio-logging" % "2.1.9",
"org.scalanlp" %% "breeze" % "2.1.0",
"org.apache.spark" %% "spark-core" % "3.3.0",
"org.apache.spark" %% "spark-sql" % "3.3.0", 
"com.github.scopt" %% "scopt" % "4.1.0"
//"org.monarchinitiative.phenol" % "phenol-core" % "2.0.0-RC4",
//"org.monarchinitiative.phenol" % "phenol-io" % "2.0.0-RC4",
//"org.monarchinitiative.phenol" % "phenol-annotations" % "2.0.0-RC4",
// "org.typelevel" %% "cats-core" % "2.9.0"
)

// mainClass in assembly := Some("com.domain.Main")

// ... Spark 
assembly / assemblyMergeStrategy := {
  case PathList("org","aopalliance", xs @ _*) => MergeStrategy.last
  case PathList("javax", "inject", xs @ _*) => MergeStrategy.last
  case PathList("javax", "servlet", xs @ _*) => MergeStrategy.last
  case PathList("javax", "activation", xs @ _*) => MergeStrategy.last
  case PathList("javax", "validation", xs @ _*) => MergeStrategy.last
  case PathList("javax", "annotation", xs @ _*) => MergeStrategy.last
  case PathList("org", "apache", xs @ _*) => MergeStrategy.last
  case PathList("com", "google", xs @ _*) => MergeStrategy.last
  case PathList("com", "esotericsoftware", xs @ _*) => MergeStrategy.last
  case PathList("com", "codahale", xs @ _*) => MergeStrategy.last
  case PathList("com", "yammer", xs @ _*) => MergeStrategy.last
  case PathList("com", "fasterxml", xs @ _*) => MergeStrategy.first
  case PathList("org", "slf4j", xs @ _*) => MergeStrategy.last
  case PathList("org", "antlr", xs @ _*) => MergeStrategy.last
  case PathList("google", "protobuf", xs @ _*) => MergeStrategy.last
  case "module-info.class" => MergeStrategy.discard
  case "git.properties" => MergeStrategy.discard
  case "about.html" => MergeStrategy.rename
  case "META-INF/ECLIPSEF.RSA" => MergeStrategy.last
  case "META-INF/mailcap" => MergeStrategy.last
  case "META-INF/mimetypes.default" => MergeStrategy.last
  case "META-INF/io.netty.versions.properties" => MergeStrategy.last  
  case "META-INF/versions/9/module-info.class" => MergeStrategy.last
  case "META-INF/org/apache/logging/log4j/core/config/plugins/Log4j2Plugins.dat" => MergeStrategy.last
  case "plugin.properties" => MergeStrategy.last
  case "log4j.properties" => MergeStrategy.last
  case "overview.html" => MergeStrategy.last  // Added this for 2.1.0 I think
  case x =>
    val oldStrategy = (assembly / assemblyMergeStrategy).value
    oldStrategy(x)
    // MergeStrategy.last
}
