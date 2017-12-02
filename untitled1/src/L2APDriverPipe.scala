import java.io.{File, FileOutputStream, PrintWriter}

import org.apache.commons.io.FileUtils
import org.apache.hadoop.hive.ql.exec.spark.session.SparkSession
import org.apache.spark._
import org.apache.spark.SparkContext._
import org.apache.spark.streaming._
import org.apache.spark.streaming.StreamingContext._

import scala.io.Source

object L2APDriverPipe {

  /** Makes sure only ERROR messages get logged to avoid log spam. */
  def setupLogging() = {
    import org.apache.log4j.{Level, Logger}
    val rootLogger = Logger.getRootLogger()
    rootLogger.setLevel(Level.ERROR)
  }

  /** Our main function where the action happens */
  def main(args: Array[String]) {



    // Set up a Spark streaming context  that runs locally using
    // all CPU cores and one-second batches of data
    val sc = new SparkContext("local[*]", "tempPipe")

    // Get rid of log spam (should be called after the context is set up)
    setupLogging()

    val thresholdV = 0.7

    val writerHeader = new PrintWriter(new FileOutputStream("/home/shalin/Downloads/pl2ap/cluHeaders.txt", false))

    writerHeader.write(thresholdV.toString)

    writerHeader.close()

    val partitions = 10

    val filename = "/home/shalin/Downloads/l2ap/build/example.clu";

    val scriptPath = "/home/shalin/Downloads/pl2ap/build/pl2ap";

    var i = 0
    var p = 1
    val writerFile = new PrintWriter(new FileOutputStream("/home/shalin/Downloads/untitled1/finalexample.csr", false))

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

    val storage = "/home/shalin/Downloads/untitled1/texts"

    val dataRDD = sc.textFile("/home/shalin/Downloads/untitled1/finalexample.csr", partitions)

    dataRDD.saveAsTextFile(storage)

    var myList = List[String]()
    if(partitions > 2){
      val k = 0;
      val j = 0;
      for(k <- 0 to partitions-1){
        for (j <- k+1 to partitions-1){
            myList = "/home/shalin/Downloads/untitled1/texts/part-0000"+k+" /home/shalin/Downloads/untitled1/texts/part-0000"+j :: myList
        }
      }
    }else{
      myList = "/home/shalin/Downloads/untitled1/texts/part-00000 /home/shalin/Downloads/untitled1/texts/part-00001" :: myList
    }

    var listpartition = 0
    if(partitions > 2){
      listpartition = partitions*partitions - partitions;
    }else{
      listpartition = 3
    }

    val finalListPartition = listpartition/2;

    val filesRDD = sc.parallelize(myList, finalListPartition)

    val pipedFiles = filesRDD.pipe(scriptPath)

    val writerOutput = new PrintWriter(new FileOutputStream("/home/shalin/Downloads/pl2ap/Output.txt", false))

    pipedFiles.collect().foreach(x => if(x.charAt(0)>='0' && x.charAt(0)<='9') writerOutput.write(x+"\n"))

    writerOutput.close()

    FileUtils.deleteDirectory(new File("/home/shalin/Downloads/untitled1/texts"))

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