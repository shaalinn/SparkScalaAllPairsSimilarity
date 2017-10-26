

import org.apache.hadoop.hive.ql.exec.spark.session.SparkSession
import org.apache.spark._
import org.apache.spark.SparkContext._
import org.apache.spark.streaming._
import org.apache.spark.streaming.twitter._
import org.apache.spark.streaming.StreamingContext._

object tempPipe {

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


    val scriptPath = "/Users/samin/Downloads/l2ap/build/apss"



    val dataRDD = sc.textFile("/Users/samin/Downloads/doc2mat-1.0/ds1.raw")

    val pipeRDD = dataRDD.pipe(scriptPath)

    pipeRDD.collect().map(println(_))
    //dataRDD.collect().map(println(_))

  }
}
