name := "clusterSpark1"

version := "0.1"

scalaVersion := "2.12.2"

resolvers ++= Seq(
  "apache-snapshots" at "http://repository.apache.org/snapshots/",
  Resolver.sonatypeRepo("public")
)

libraryDependencies += "org.apache.spark" % "spark-core_2.11" % "2.2.0" % "provided"

libraryDependencies += "commons-io" % "commons-io" % "2.4"

assemblyJarName in assembly := s"sparkcluster.jar"

mainClass in assembly := Some("com.cmpe295b.Main")





assemblyOption in assembly :=
  (assemblyOption in assembly).value.copy(includeScala=false)

