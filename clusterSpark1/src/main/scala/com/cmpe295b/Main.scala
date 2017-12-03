package com.cmpe295b
import java.io.{File, FileOutputStream, PrintWriter}

import org.apache.commons.io.FileUtils
import org.apache.spark._

import scala.io.Source

object Main {
  /** Makes sure only ERROR messages get logged to avoid log spam. */
  def setupLogging() = {
    import org.apache.log4j.{Level, Logger}
    val rootLogger = Logger.getRootLogger()
    rootLogger.setLevel(Level.ERROR)
  }

  /** Our main function where the action happens */
  def main(args: Array[String]) {

    val BINARY_COMMON_LOCATION = "/home/hadoop/SparkScalaAllPairsSimilarity"
    val STORAGE_COMMON_LOCATION = "/home/hadoop/efs/storage_temp"
    val OUTPUT_COMMON_LOCATION = "/home/hadoop/efs/storage_temp"


    // Set up a Spark streaming context  that runs locally using
    // all CPU cores and one-second batches of data

    val conf1 = new SparkConf()
    conf1.set("spark.hadoop.validateOutputSpecs", "false")
    conf1.set("fs.default.name","file:///")
    conf1.setAppName("tempPipe")
    conf1.setMaster("local[*]")


    val sc2 = new SparkContext(conf1)



    // Get rid of log spam (should be called after the context is set up)
    setupLogging()

    val thresholdV = args{1}.toInt

    val writerHeader = new PrintWriter(new FileOutputStream(STORAGE_COMMON_LOCATION+"/cluHeaders.txt", false))

    writerHeader.write(thresholdV.toString)

    writerHeader.close()

    val partitions = args{0}.toFloat

    val filename = STORAGE_COMMON_LOCATION+"/"+args{2};

    val scriptPath = BINARY_COMMON_LOCATION+"/pl2ap/build/pl2ap";

    var i = 0
    var p = 1
    val writerFile = new PrintWriter(new FileOutputStream(STORAGE_COMMON_LOCATION+"/finalexample.csr", false))

    if(filename.substring(filename.length-3,filename.length).equals("clu")){
      for(line <- Source.fromFile(filename).getLines()){
        i = i + 1;
        if(line.split(" ").length ==3 && i == 1){

        }else{
          writerFile.write(i+line+"\n")
        }
      }
    }else if(filename.substring(filename.length-3,filename.length).equals("csr")){
      for(line <- Source.fromFile(filename).getLines()){
        p = p + 1;
        writerFile.write(i+line+"\n")
      }
    }

    val storage = STORAGE_COMMON_LOCATION+"/texts"

    val dataRDD = sc2.textFile("file:///"+STORAGE_COMMON_LOCATION+"/finalexample.csr", partitions)

    dataRDD.saveAsTextFile("file:///"+storage)

    sc2.stop()

    val conf = new SparkConf()
    conf.set("spark.hadoop.validateOutputSpecs", "false")
    conf.set("fs.default.name","file:///")
    conf.setAppName("tempPipe")

    val sc = new SparkContext(conf)

    var myList = List[String]()
    if(partitions > 2){
      var k = 0
      while(k < partitions){
        var j = k+1
        while (j < partitions){
          myList = STORAGE_COMMON_LOCATION+"/texts/part-0000"+k+" "+STORAGE_COMMON_LOCATION+"/texts/part-0000"+j :: myList
          j = j + 1
        }
        k = k + 1
      }
    }else{
      myList = STORAGE_COMMON_LOCATION+"/texts/part-00000 "+STORAGE_COMMON_LOCATION+"/texts/part-00001" :: myList
    }

    var listpartition = 0
    if(partitions > 2){
      listpartition = partitions*partitions - partitions;
    }else{
      listpartition = 3
    }

    //println(myList)

    val finalListPartition = listpartition / 2;

    //println(finalListPartition)

    val filesRDD = sc.parallelize(myList, finalListPartition)

    //filesRDD.saveAsTextFile("file:///"+STORAGE_COMMON_LOCATION+"/texts")

    val pipedFiles = filesRDD.pipe(scriptPath)

    val writerOutput = new PrintWriter(new FileOutputStream(OUTPUT_COMMON_LOCATION+"/Output.txt", false))

    pipedFiles.collect.foreach( x => if(x.charAt(0)>='0' && x.charAt(0)<='9') writerOutput.write(x+"\n"))


    //pipedFiles.collect().foreach(x => if(x.charAt(0)>='0' && x.charAt(0)<='9') writerOutput.write(x+"\n"))

    writerOutput.close()

    //FileUtils.deleteDirectory(new File(STORAGE_COMMON_LOCATION+"/texts"))

    /*val originalRDD = dataRDD.cache()

    val header = dataRDD.first()

    val xyz = header.split(" ")

    val ogData = originalRDD.filter(row => row != header)

    println(xyz(0)+ xyz(1))

    val writer = new PrintWriter(new FileOutputStream("/home/shalin/Downloads/l2ap/build/cluHeaders.txt", false))

    writer.write(xyz(1)+" "+xyz(2)+"\n"+thresholdV)

    writer.close()

    val pipeRDD = ogData.pipe(scriptPath)

    val writer1 = new PrintWriter(new FileOutputStream("/home/shalin/Desktop/finalOutput.txt", false))

    pipeRDD.collect().foreach(x => writer1.write(x+"\n"))

    writer1.close()*/

    //dataRDD.collect().map(println(_))

  }
}
