package com.cmpe295b
import java.io.{File, FileOutputStream, PrintWriter}

import org.apache.commons.io.FileUtils
import org.apache.spark._

import scala.io.Source
import scala.util.Try

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

    val thresholdV = args{1}.toFloat

    val writerHeader = new PrintWriter(new FileOutputStream(STORAGE_COMMON_LOCATION+"/cluHeaders.txt", false))

    writerHeader.write(thresholdV.toString)

    writerHeader.close()

    //TODO dynamic partitions in next two lines
    val partitions = args{0}.toInt

    println("partitions are "+partitions)

    val filename = STORAGE_COMMON_LOCATION+"/"+args{2};

    val scriptPath = BINARY_COMMON_LOCATION+"/pl2ap/build/pl2ap";
    val scriptPathSingle = BINARY_COMMON_LOCATION+"/pl2apPL2AP/build/pl2ap";

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

    //val dataRDD = dataRDD1.repartition(partitions)

    println(dataRDD.partitions.length)

    val realParitions = dataRDD.partitions.length

    dataRDD.saveAsTextFile("file:///"+storage)

    sc2.stop()

    val conf = new SparkConf()
    conf.set("spark.hadoop.validateOutputSpecs", "false")
    conf.set("fs.default.name","file:///")
    conf.setAppName("tempPipe")


    val sc = new SparkContext(conf)

    var k = 0;
    var j = 0;
    var myList = List[String]()
    if(realParitions > 2){
      while(k < realParitions){
        j = k+1
        while (j < realParitions){
          myList = STORAGE_COMMON_LOCATION+"/texts/part-"+f"$k%05d"+" "+STORAGE_COMMON_LOCATION+"/texts/part-"+f"$j%05d" :: myList
          j = j + 1
        }
        k = k + 1
      }
    }else{
      myList = STORAGE_COMMON_LOCATION+"/texts/part-00000 "+STORAGE_COMMON_LOCATION+"/texts/part-00001" :: myList
    }

    println("size of myList is "+ myList.length)

    var listpartition = 0
    if(realParitions > 2){
      listpartition = realParitions*realParitions - realParitions;
    }else{
      listpartition = 3
    }

    if(realParitions>1){
      var a = 0;
      while (a < realParitions){
        myList = STORAGE_COMMON_LOCATION+"/texts/part-"+f"$a%05d"+" "+STORAGE_COMMON_LOCATION+"/texts/part-"+f"$a%05d" :: myList
        a = a + 1
      }
    }
    
    println("my list "+myList)

    println("final list partitions "+listpartition)

    val finalListPartition = listpartition / 2;

    val finalListPartition1 = finalListPartition + realParitions

    println("final list partitions "+finalListPartition)

    val filesRDD = sc.parallelize(myList, finalListPartition1)

    println("partition made for list "+ filesRDD.partitions.length)

    //filesRDD.saveAsTextFile("s3n://sparkcluster-cmpe295b/")
    //filesRDD.saveAsTextFile("file:///"+STORAGE_COMMON_LOCATION+"/texts")

    val pipedFiles = filesRDD.pipe(scriptPath)


    //fixme single file pl2ap
    var q = 0
    var myListSingle = List[String]()
    while(q < realParitions-1){
      myListSingle = STORAGE_COMMON_LOCATION+"/texts/part-"+f"$q%05d" :: myListSingle
      q = q + 1
    }

    println("combination "+ q)
    println("size of myListSingle is "+ myListSingle.length)

    println("myListSingle is "+ myListSingle)

    val singleFilesRDD = sc.parallelize(myListSingle, realParitions)

    val singleFilesRDD1 = singleFilesRDD.repartition(realParitions)

    val singlePipedFiles = singleFilesRDD1.pipe(scriptPathSingle)

    //val finalResultRDD = pipedFiles.union(singlePipedFiles)

    //val writerOutput = new PrintWriter(new FileOutputStream(OUTPUT_COMMON_LOCATION+"/Output.txt", false))

    //pipedFiles.collect.foreach( x => if(x.charAt(0)>='0' && x.charAt(0)<='9') writerOutput.write(x+"\n"))

    //pipedFiles.coalesce(1, true).saveAsTextFile("s3n://sparkcluster-cmpe295b/")

    pipedFiles.saveAsTextFile("s3n://sparkcluster-cmpe295b/")
     //finalResultRDD.coalesce(realParitions, true).saveAsTextFile("s3n://sparkcluster-cmpe295b/")




    //pipedFiles.collect().foreach(x => if(x.charAt(0)>='0' && x.charAt(0)<='9') writerOutput.write(x+"\n"))

    //writerOutput.close()

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
