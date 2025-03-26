
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigInteger;
import java.security.InvalidAlgorithmParameterException;
import java.security.InvalidKeyException;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.Random;
import java.util.Vector;
import java.util.concurrent.*;
import java.lang.management.ManagementFactory;
import com.sun.management.OperatingSystemMXBean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;


public class Asymmetric_DH_together{
	public static void main(String[] args) {
		try {
			Experiment.measureTime();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}

class Experiment{
	public static void measureTime() throws IOException {
		//String OUTPUT = "ADH_k=24SEDH_2048_d4.csv";
		int NUM_RUN = 110;
		int NUM_CONDITION = 1;
		//int bitSize_m;
		int keysize = 16384;  // 既存のkeysize
		Random random = new SecureRandom();
		BigInteger p = BigInteger.probablePrime(keysize, random);
		BigInteger g = randomBInt(random, p);

		Scanner scanner = new Scanner(System.in);  // 標準入力をスキャナで読み込む
		PrintWriter pw = null;  // ファイル出力用のPrintWriter

        // 2. 複数のdをループで処理する
    while (scanner.hasNextInt()) {  // 標準入力に次の整数がある限り実行
        int d = scanner.nextInt();  // dを1つ読み込む
        System.out.println("d = " + d);  // 現在のdを表示


        int bitSize_m = keysize / d;  // mの計算 (keysize / d)
        // keysize と d を文字列に反映
				String OUTPUT = String.format("WithBobADH_k=24SEDH_%d_bit_d_%d_m_%d.csv", keysize, d, bitSize_m);
			  pw = new PrintWriter(new BufferedWriter(new FileWriter(OUTPUT, true)));  // 追記モードでファイルを開く

			 // 4. 実際の処理 (ここに実験内容を追加します)
			  System.out.printf("現在の設定: keysize = %d, d = %d, m = %d\n", keysize, d, bitSize_m);

	//	PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(OUTPUT)));
				Experiment.printHeader(pw);



				double averageCpuLoad =0;
				double coresUsed = 0;

				for (int i = 0; i < NUM_CONDITION; i++) {
						//int keysize = 2048;
						//int d = 6;
						//bitSize_m = keysize/d;
						System.out.println("size:" + keysize);
						System.out.println("m=" + bitSize_m);

						BigInteger m = BigInteger.valueOf(bitSize_m);

						System.out.println("d=" + d);



				//		BigInteger[][] yB1List = new BigInteger[NUM_RUN][d];
		    //		BigInteger[][] yB2List = new BigInteger[NUM_RUN][d];
				//		BigInteger[][] baseList = new BigInteger[NUM_RUN][d];
				//		BigInteger[] bList = new BigInteger[NUM_RUN]; // xBも事前に生成して保持


					//	System.out.printf("Bobの公開鍵の事前計算中\n");
						// CPU使用率を取得するための設定
						OperatingSystemMXBean osBean = ManagementFactory.getPlatformMXBean(OperatingSystemMXBean.class);
						int availableProcessors = Runtime.getRuntime().availableProcessors();
						List<Double> cpuLoadSamples = new ArrayList<>();


				    // 事前に NUM_RUN 分の yB1 と yB2 を生成
				   // for (int j = 0; j < NUM_RUN; j++) {

				     //   BigInteger[] base = createBaseVec(m, d, p.subtract(BigInteger.ONE), BigInteger.ONE);
						//		BigInteger b =  randomBInt(random, p);
						//		BigInteger[] baseB = createBaseVec(m, d, p.subtract(BigInteger.ONE), b);

						//		baseList[j] = base;
						//		bList[j] = b;
		        //		yB1List[j] = schurExpVec(baseB, g, p, d);
		        //		yB2List[j] = schurExpVec(base, g, p, d); // NB[0]を使ってyB2を生成
				    //}
						//System.out.printf("Bobの公開鍵の事前計算終了\n");


						for (int j = 0; j < NUM_RUN; j++) {
								// CPU使用率サンプリングスレッド
								Thread cpuMonitor = new Thread(() -> {
										try {
												while (!Thread.currentThread().isInterrupted()) {
														double processCpuLoad = osBean.getProcessCpuLoad();
														if (processCpuLoad >= 0) {
																cpuLoadSamples.add(processCpuLoad);
														}
														Thread.sleep(500);  // 500ミリ秒ごとにサンプリング
												}
										} catch (InterruptedException e) {
												Thread.currentThread().interrupt();
										}
								});

								// CPUモニタースレッド開始
								cpuMonitor.start();
								long begin = System.currentTimeMillis();

								try {
										// 実際の処理（ここに対象の処理を記述）

										printResurt(pw, d, p, bitSize_m, run(d, p, g, bitSize_m),averageCpuLoad,coresUsed);
								} catch (InvalidKeyException | NoSuchAlgorithmException | InvalidAlgorithmParameterException e) {
										e.printStackTrace();
								}

								long end = System.currentTimeMillis();

								// CPUモニタースレッドを停止
								cpuMonitor.interrupt();
								try {
										cpuMonitor.join();  // CPUモニタースレッドの終了を待つ
								} catch (InterruptedException e) {
										Thread.currentThread().interrupt();
								}

								// 平均CPU使用率を計算
								double totalCpuLoad = 0;
								for (double load : cpuLoadSamples) {
										totalCpuLoad += load;
								}
								averageCpuLoad = cpuLoadSamples.isEmpty() ? 0 : totalCpuLoad / cpuLoadSamples.size();
								coresUsed = averageCpuLoad * availableProcessors;
								// CPU使用率の出力
				        System.out.printf("平均CPU使用率: %.2f%%, 使用されたコア数: %.2f/%d%n", averageCpuLoad * 100, coresUsed, availableProcessors);

						}

				}

				pw.flush();
		}
			// 5. ファイルとスキャナを閉じる
	  if (pw != null) {
	     pw.close();  // ファイル出力用のPrintWriterを閉じる
	  }
	  scanner.close();  // Scannerを閉じる
	  System.out.println("すべてのdの実行が終了しました");
	}


	private static long[] run(int d, BigInteger p, BigInteger g, int bitSize_m)
			throws NoSuchAlgorithmException, InvalidAlgorithmParameterException, InvalidKeyException {

		Random random = new SecureRandom();
		BigInteger[][] data;
		long[] time = new long[6];

		// Master Key pair
		//BigInteger q = p.subtract(BigInteger.ONE);
		//BigInteger[][] xB = createMatrix(random, p, d, p);
		//BigInteger[][][] NB = makeInvertibleMatrix(random, d, q);
		BigInteger m = BigInteger.valueOf(bitSize_m);

		BigInteger[] base = createBaseVec(m, d, p.subtract(BigInteger.ONE), BigInteger.ONE);
		BigInteger b =  randomBInt(random, p);
		BigInteger[] baseB = createBaseVec(m, d, p.subtract(BigInteger.ONE), b);

		int int_two = 2;
		BigInteger two = BigInteger.valueOf(int_two);
		BigInteger[] xA = createVec(random, two.pow(bitSize_m), d, p);

		time[0] = System.nanoTime();

		BigInteger[] yB1 = schurExpVec(baseB, g, p, d);
		BigInteger[] yB2 = schurExpVec(base, g, p, d);
		//BigInteger[][] yB1 = schurExp(xB, g, p, d);
		//BigInteger[][] yB2 = schurExp(NB[0], g, p, d);

		time[1] = System.nanoTime();

		ompSEVec pse = new ompSEVec(yB2, xA, p, d);

		pse.ompSERightVec();
		BigInteger yA = pse.getResult();
		//System.out.println("calc_finished_yA" + "\n");
		time[2] = System.nanoTime();

		ompSEVec pse2 = new ompSEVec(yB1, xA, p, d);

		pse2.ompSERightVec();
		BigInteger sssk = pse2.getResult();
		//System.out.println("calc_finished_KA" + "\n");

		//System.out.println("k_A=");
	//	System.out.println(sssk.toString());
		//sssk.print();


		// Master SSK
		time[3] = System.nanoTime();

	  //BigInteger[][] msk = multi(xB, NB[1],  p.subtract(BigInteger.ONE), d);

	  //ompSEVecLeft pse3 = new ompSEVecLeft(yA, msk, p, d);
	//	pse3.ompSELeftVec();
 		BigInteger mssk = yA.modPow(b, p);

	//	System.out.println("k_B=");
	//	System.out.println(mssk.toString());
//		mssk.print();


//		PararllelSE pse3_test = new PararllelSE(yA, msk, p, d,64);
//		pse3_test.SELeft();
//		ModularMatrix mssk_test = pse3_test.getResult();
//		System.out.println("k_B_test=");
//		mssk_test.print();

		time[4] = System.nanoTime();

		for (int i = 0; i < 5; i++) {
        time[i] = time[i] / 1_000; // ns → μs
    }

		// **スレッドごとの平均時間を計算**
    List<Long> executionTimes = pse.getExecutionTimes();
    double avgThreadExecutionTime = executionTimes.stream().mapToLong(Long::longValue).average().orElse(0) / 1_000.0; // マイクロ秒換算
    time[5] = (long) avgThreadExecutionTime; // avgThreadExecutionTime を追加

		System.out.printf("合計実行時間: %.2f μs, 平均スレッド実行時間: %.2f μs\n",
                  (double) (time[4] - time[0]), avgThreadExecutionTime);
    // **結果を出力**
    System.out.printf("平均スレッド実行時間: %.2f μs\n", avgThreadExecutionTime);


//			for (int j = 0; j < d; j++) {
		if (mssk.subtract(sssk) != BigInteger.ZERO){
			throw new RuntimeException("Keys are not equals.");
		}
//		}

		return time;
	}


	private static void printHeader(PrintWriter pw) {
		String[] label = { "dim", "bit.length.p", "bit.length.m", "MKP", "SKP", "SSK", "MSK",
                       "TotalMaster", "TotalSlave", "averageCpuLoad", "coresUsed", "AvgThreadExecTime"};
    String str = String.join(", ", label) + "\n";
    pw.print(str);
	}

	private static void printResurt(PrintWriter pw, int d, BigInteger p, int bitSize_m,  long[] time, double averageCpuLoad, double coresUsed) {
		double avgThreadExecutionTime = time[5]; // run() から取得

    String str = d + ", " + p.bitLength() + ", " + bitSize_m + ", ";
    str += (time[1] - time[0]) + ", "; // MKP
    str += (time[2] - time[1]) + ", "; // SKP
    str += (time[3] - time[2]) + ", "; // SSK
    str += (time[4] - time[3]) + ", "; // MSK
    str += (time[1] - time[0] + time[4] - time[3]) + ", "; // MasterTotal
    str += (time[3] - time[1]) + ", "; // SlaveTotal
    str += averageCpuLoad + ", ";
    str += coresUsed + ", ";
    str += avgThreadExecutionTime + "\n"; // avgThreadExecutionTime を追加

    pw.print(str);
		System.out.println("[" + p.bitLength() + " bits]:" + (time[2] - time[1]) + " μs for Alice y_A,");
    System.out.println("[" + p.bitLength() + " bits]:" + (time[3] - time[2]) + " μs for Alice _A,");
	}


	private static int condD(int cond, int swt) {
		if ((swt & 0x100) > 0)
			return (cond + 2);
		else
			return 10;
	}

	private static BigInteger condP(int cond, int swt) {
		if ((swt & 0x010) > 0)
			return BigInteger.probablePrime((cond + 1) * 8 - 1, new Random());
		else
			return new BigInteger("2147483647");
	}

	private static int condI(int cond, int swt) {
		if ((swt & 0x001) > 0)
			return 20 * (cond + 1);
		else
			return 10;
	}

	public static String toString(BigInteger[][] a, int d) {
		String res = "\r\n";
		String s = "[";
		for (int i = 0; i < d; i++) {
			for (int j = 0; j < d; j++) {
				s = s.concat(a[i][j].toString());
				if (j < d - 1)
					s += ", ";
				else if (i == d - 1)
					s += "]";
			}
			res += s + "\r\n";
			s = " ";
		}
		return res;
	}

	public static String toStringVec(BigInteger[] a, int d) {
		String res = "\r\n";
		String s = "[";
		for (int i = 0; i < d; i++) {
			s = s.concat(a[i].toString());
			if (i < d - 1)
				s += ", ";
			else if (i == d - 1)
				s += "]";

			res += s + "\r\n";
			s = " ";
		}
		return res;
	}

	public static BigInteger randomBInt(Random random, BigInteger modulus) {
		if (modulus.compareTo(new BigInteger("2")) < 0)
			throw new IllegalArgumentException("mod must be gleater or equalsto 2.");

		BigInteger res;
		BigInteger MAX = new BigInteger(Integer.toString(Integer.MAX_VALUE));

		if (modulus.compareTo(MAX) > 0) {
			BigInteger m = modulus;
			Vector<Integer> v = new Vector<Integer>();

			// Split biginteger into size of integer.
			while (!m.equals(BigInteger.ZERO)) {
				v.add(random.nextInt(m.mod(MAX).intValue()));
				m = m.divide(MAX);
			}

			res = BigInteger.ZERO;
			int cnt = 0;

			// Assembling the cut value again
			for (Integer b : v) {
				res = res.add(new BigInteger(b.toString()).multiply(MAX.pow(cnt)));
				cnt++;
			}

			return res;
		} else if (modulus.compareTo(BigInteger.ZERO) > 0) {
			// In case of mod < IntMAX
			do
				res = new BigInteger(Integer.toString(random.nextInt(modulus.intValue())));
			while (res.equals(BigInteger.ZERO));
			return res;
		} else
			throw new IllegalArgumentException("'mod' must be greater than 0.");
	}

	public static BigInteger randomInvBInt(Random random, BigInteger modulus) {
		if (modulus.compareTo(new BigInteger("2")) < 0)
			throw new IllegalArgumentException("mod must be gleater or equalsto 2.");
		while (true) {
			try {
				return randomBInt(random, modulus).modInverse(modulus);
			} catch (ArithmeticException e) {
				// not invertible. try again.
			}
		}
	}

	public static BigInteger[][] diag(BigInteger val, int dim) {
		if (dim <= 0)
			throw new IllegalArgumentException("\'dim\' must be positive.(" + dim + ")");
		BigInteger[][] res = new BigInteger[dim][dim];

		for (int i = 0; i < dim; i++)
			for (int j = 0; j < dim; j++)
				res[i][j] = (i == j) ? val : BigInteger.ZERO;
		return res;
	}

	public static BigInteger[][] idm(int dim) {
		return diag(BigInteger.ONE, dim);
	}

	public static BigInteger[][] createMatrix(Random random, BigInteger m, int dim, BigInteger mod) {
		BigInteger[][] res = new BigInteger[dim][dim];

		for (int i = 0; i < dim; i++)
			for (int j = 0; j < dim; j++)
				res[i][j] = randomBInt(random, m);

		return res;
	}

	public static BigInteger[] createVec(Random random, BigInteger m, int dim, BigInteger mod) {
		BigInteger[] res = new BigInteger[dim];

		for (int i = 0; i < dim; i++)
			res[i] = randomBInt(random, m);

		return res;
	}
//new 202502
	public static BigInteger[] createBaseVec(BigInteger m, int dim, BigInteger mod, BigInteger b) {
		BigInteger[] res = new BigInteger[dim];
		BigInteger two = BigInteger.TWO;
		BigInteger exponent;
		BigInteger t;
		BigInteger tmp;
		for (int i = 0; i < dim; i++){
			t = BigInteger.valueOf(i);
			exponent = m.multiply(t);
			tmp = two.modPow(exponent, mod);
			res[i] = tmp.multiply(b).mod(mod);
			//System.out.println("test: " + res[i]  + "\n");
		}

		return res;

	}

	public static BigInteger[][] createMatrixInvertibleLTM(Random random, int dim, BigInteger mod) {
		BigInteger[][] res = new BigInteger[dim][dim];

		for (int i = 0; i < dim; i++)
			for (int j = 0; j < dim; j++)
				if (i == j)
					res[i][j] = randomInvBInt(random, mod);
				else if (i < j)
					res[i][j] = BigInteger.ZERO;
				else
					res[i][j] = randomBInt(random, mod);

		return res;
	}

	public static BigInteger[][] createMatrixInvertibleUTM(Random random, int dim, BigInteger mod) {
		BigInteger[][] res = new BigInteger[dim][dim];

		for (int i = 0; i < dim; i++)
			for (int j = 0; j < dim; j++)
				if (i == j)
					res[i][j] = randomInvBInt(random, mod);
				else if (i > j)
					res[i][j] = BigInteger.ZERO;
				else
					res[i][j] = randomBInt(random, mod);

		return res;
	}

	public static BigInteger[][][] makeInvertibleMatrix(Random random, int dim, BigInteger modulus) {
		if (modulus.compareTo(new BigInteger("2")) < 0)
			throw new IllegalArgumentException("mod must be gleater or equalsto 2.");

		BigInteger[][] l = createMatrixInvertibleLTM(random, dim, modulus);
		BigInteger[][] u = createMatrixInvertibleUTM(random, dim, modulus);

//		System.out.println("u=");
//		System.out.println(toString(u,dim));
//
//		System.out.println("u^t=");
//		System.out.println(toString(transpose(u,dim),dim));

		// generate LTM inverse
		BigInteger[][] li = generateLTMInv(l, modulus, dim);

		// generate UTM inverse
		BigInteger[][] ui = transpose(generateLTMInv(transpose(u,dim), modulus,dim),dim);

		BigInteger[][][] res = { multi(l,u,modulus, dim), multi(ui,li,modulus,dim) };
		return res;
	}

	private static BigInteger[][] generateLTMInv(BigInteger[][] ltm, BigInteger mod, int dim ) {

		BigInteger[][] linv = new BigInteger[dim][dim];

		// set 1st row
		for (int j = 0; j < dim; j++) {
			if (0 == j)
				linv[0][j] = BigInteger.ONE;
			else
				linv[0][j] = BigInteger.ZERO;
		}

		for (int i = 0; i < dim; i++)
			for (int j = 0; j < dim; j++) {
				// divide ith row of L' by L[i,i]
				linv[i][j] = linv[i][j].multiply(ltm[i][i].modInverse(mod));

				if (i < dim - 1) {
					if (i + 1 == j)
						linv[i + 1][j] = BigInteger.ONE;
					else
						linv[i + 1][j] = BigInteger.ZERO;

					for (int k = j; k < i + 1; k++)
						// L'[i+1,j] = L'[i+1,j] - L[i+1,k] * L'[k,j]
						linv[i + 1][j] = linv[i + 1][j].subtract(ltm[i+1][k].multiply(linv[k][j]));
					//	linv[i + 1][j] = linv[i + 1][j].mod(mod) ;
				}
			}
		return linv;
	}

	public static BigInteger[][] transpose(BigInteger[][] a, int dim) {
		BigInteger[][] res = new BigInteger[dim][dim];
		for (int i = 0; i < dim; i++)
			for (int j = 0; j < dim; j++)
				res[i][j] = a[j][i];
		return res;
	}

	public static BigInteger[][] multi(BigInteger[][] a, BigInteger[][] b, BigInteger mod, int d) {
		BigInteger[][] res = new BigInteger[d][d];
		for (int i = 0; i < d; i++)
			for (int j = 0; j < d; j++) {
				BigInteger v = BigInteger.ZERO;
				for (int k = 0; k < d; k++)
					v = v.add(a[i][k].multiply(b[k][j]));
				res[i][j] = v.mod(mod);
			}
		return res;
	}

	public static BigInteger schurExp(BigInteger base, BigInteger exponent, BigInteger modulus) {
		if (modulus.compareTo(new BigInteger("2")) < 0)
			throw new IllegalArgumentException("mod must be gleater or equalsto 2.");

		if (base.equals(BigInteger.ZERO))
			return BigInteger.ZERO;
		else
			return base.modPow(exponent, modulus);
	}

	public static BigInteger[][] schurExp(BigInteger[][] mat, BigInteger base, BigInteger mod, int dim) {
        BigInteger[][] res = new BigInteger[dim][dim];



        // 並列化のためのExecutorServiceを使用
				int availableProcessors = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(availableProcessors);
        List<Future<?>> futures = new ArrayList<>();

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                int row = i;
                int col = j;
                Future<?> future = executor.submit(() -> {
                    res[row][col] = schurExp(base, mat[row][col], mod);
                });
                futures.add(future);
            }
        }

        // 全てのタスクが完了するまで待機
        for (Future<?> future : futures) {
            try {
                future.get();
            } catch (InterruptedException | ExecutionException e) {
                e.printStackTrace();
            }
        }

        // ExecutorServiceを終了
        executor.shutdown();


        return res;
    }


		public static BigInteger[] schurExpVec(BigInteger[] vec, BigInteger base, BigInteger mod, int dim) {
	        BigInteger[] res = new BigInteger[dim];



	        // 並列化のためのExecutorServiceを使用
					int availableProcessors = Runtime.getRuntime().availableProcessors();
	        ExecutorService executor = Executors.newFixedThreadPool(availableProcessors);
	        List<Future<?>> futures = new ArrayList<>();

	        for (int i = 0; i < dim; i++) {

	            int row = i;

	            Future<?> future = executor.submit(() -> {
	                res[row] = schurExp(base, vec[row], mod);
	            });
	            futures.add(future);
	        }

	        // 全てのタスクが完了するまで待機
	        for (Future<?> future : futures) {
	            try {
	                future.get();
	            } catch (InterruptedException | ExecutionException e) {
	                e.printStackTrace();
	            }
	        }

	        // ExecutorServiceを終了
	        executor.shutdown();


	        return res;
	    }

}

class ompSEVec{
	BigInteger[] a;
	BigInteger[] b;
	BigInteger c;
	BigInteger[] tmp;
	BigInteger data;
	BigInteger p;
	int MATRIX_SIZE;
	private List<Long> executionTimes = new ArrayList<>(); // 実行時間を保存するリスト

	ompSEVec(BigInteger[] a, BigInteger[] b, BigInteger p, int d) {
		// assumption : a and b are both double[MATRIX_SIZE][MATRIX_SIZE]
		this.a = a;
		this.b = b;
		this.MATRIX_SIZE = d;
		this.p = p;
		this.data = BigInteger.ONE;
		this.tmp = new BigInteger[MATRIX_SIZE];


		for (int i = 0; i < MATRIX_SIZE; ++i) {
				tmp[i] = BigInteger.ONE;
		}
	}

	public void ompSERightVec() {
	        // CPU使用率を取得するための設定
	       // OperatingSystemMXBean osBean = ManagementFactory.getPlatformMXBean(OperatingSystemMXBean.class);
	        int availableProcessors = Runtime.getRuntime().availableProcessors();
	       // List<Double> cpuLoadSamples = new ArrayList<>();

	        // CPU使用率サンプリングスレッド
	     //   Thread cpuMonitor = new Thread(() -> {
	      //      try {
	        //        while (!Thread.currentThread().isInterrupted()) {
	          //          double processCpuLoad = osBean.getProcessCpuLoad();
	            //        if (processCpuLoad >= 0) {
	              //          cpuLoadSamples.add(processCpuLoad);
	                //    }
	                  //  Thread.sleep(500);  // 500ミリ秒ごとにサンプリング
	 //               }
	   //         } catch (InterruptedException e) {
	     //           Thread.currentThread().interrupt();
	       //     }
	     //   });

	        // CPUモニタースレッド開始
	      //  cpuMonitor.start();

	        // 並列化のためのExecutorServiceを使用
	        ExecutorService executor = Executors.newFixedThreadPool(availableProcessors);
	        List<Future<?>> futures = new ArrayList<>();
					executionTimes.clear();

					// 行と列の計算を並列化
	        for (int i = 0; i < MATRIX_SIZE; ++i) {

	            int row = i;

	            Future<?> future = executor.submit(() -> {
								long startTime = System.nanoTime();
                tmp[row] = modPow(a[row], b[row], p);
                long endTime = System.nanoTime();
                synchronized (executionTimes) {
                    executionTimes.add(endTime - startTime);
                }
	               // System.out.println("now calc [" + row + "]" + "\n");
	            });
	            futures.add(future);

	        }

	        // 全てのタスクが完了するまで待機
	        for (Future<?> future : futures) {
	            try {
	                future.get();
	            } catch (InterruptedException | ExecutionException e) {
	                e.printStackTrace();
	            }
	        }

	        // 行全体の計算も並列化
	        futures.clear();
					for (int i = 0; i < MATRIX_SIZE; ++i) {
					    int row = i;
					        //Future<?> future = executor.submit(() -> {
					    data = data.multiply(tmp[row]).mod(p);
					        //});
					        //futures.add(future);
					}

	        // 全ての行に対する計算が完了するまで待機
	        for (Future<?> future : futures) {
	            try {
	                future.get();
	            } catch (InterruptedException | ExecutionException e) {
	                e.printStackTrace();
	            }
	        }

	        // スレッドプールをシャットダウン
	        executor.shutdown();
	        // CPUモニタースレッドを停止
	       // cpuMonitor.interrupt();
	      //  try {
	       //     cpuMonitor.join();  // CPUモニタースレッドの終了を待つ
	      //  } catch (InterruptedException e) {
	       //     Thread.currentThread().interrupt();
	       // }

	        // 平均CPU使用率を計算
	      //  double totalCpuLoad = 0;
	      //  for (double load : cpuLoadSamples) {
	      //      totalCpuLoad += load;
	      //  }
	      //  double averageCpuLoad = cpuLoadSamples.isEmpty() ? 0 : totalCpuLoad / cpuLoadSamples.size();
	      //  double coresUsed = averageCpuLoad * availableProcessors;

	        // CPU使用率の出力
	      //  System.out.printf("平均CPU使用率: %.2f%%, 使用されたコア数: %.2f/%d%n", averageCpuLoad * 100, coresUsed, availableProcessors);
	    }


	public BigInteger getResult() {
		// for (int i = 0; i < MATRIX_SIZE; ++i) {
		// 	for (int k = 0; k < MATRIX_SIZE; ++k) {
		// 		data[i] = data[i].multiply(tmp[i][k]).mod(p);
		// 	}
		// }


		// for (int i = 0; i < MATRIX_SIZE; ++i) {
		// 	// for (int k = 0; k < MATRIX_SIZE; ++k) {
		// 	///	data[i] =  data[i].multiply(tmp[i][k]).mod(p);
		// 		data[i] =  data[i];
		// 	//}
		// }
		c = data;
		return c;
	}
	public List<Long> getExecutionTimes() {
        return executionTimes;
  }

	private static BigInteger modPow(BigInteger base, BigInteger exponent, BigInteger mod) {
		if (base.equals(BigInteger.ZERO))
			return BigInteger.ZERO;
		else
			return base.modPow(exponent, mod);
	}

}

class ompSEVecLeft{
	BigInteger[] a;
	BigInteger[][] b;
	BigInteger[] c;
	BigInteger [] data;
	BigInteger [][] tmp;
	BigInteger p;
	int MATRIX_SIZE;
	int availableProcessors = Runtime.getRuntime().availableProcessors();
	ExecutorService executor = Executors.newFixedThreadPool(availableProcessors);
	List<Future<?>> futures = new ArrayList<>();

	ompSEVecLeft(BigInteger[] a, BigInteger[][] b, BigInteger p, int d) {
		// assumption : a and b are both double[MATRIX_SIZE][MATRIX_SIZE]
		this.a = a;
		this.b = b;
		this.MATRIX_SIZE = d;
		this.p = p;
		this.data = new BigInteger[MATRIX_SIZE];
		this.tmp = new BigInteger[MATRIX_SIZE][MATRIX_SIZE];

		for (int i = 0; i < MATRIX_SIZE; ++i) {
			data[i] = BigInteger.ONE;
		}

		for (int i = 0; i < MATRIX_SIZE; ++i) {
			for (int k = 0; k < MATRIX_SIZE; ++k) {
				tmp[i][k] = BigInteger.ONE;
			}
		}
	}

	public void ompSELeftVec() {

		// 行（i）と列（k）の計算を並列化
        for (int i = 0; i < MATRIX_SIZE; ++i) {
            int row = i;
            for (int k = 0; k < MATRIX_SIZE; ++k) {
                int col = k;

                Future<?> future = executor.submit(() -> {
                    // tmp[i][k] = modPow(a[k], b[i][k], p);
                    tmp[row][col] = modPow(a[col], b[row][col], p);

                    // data[i] の更新
                    synchronized (data) {
                        data[row] = data[row].multiply(tmp[row][col]);
                    }
                });
                futures.add(future);
            }
        }

        // 全てのタスクが完了するまで待機
        for (Future<?> future : futures) {
            try {
                future.get();
            } catch (InterruptedException | ExecutionException e) {
                e.printStackTrace();
            }
        }

        // 行ごとに mod p の計算を行う（これも並列化可能）
        futures.clear();
        for (int i = 0; i < MATRIX_SIZE; ++i) {
            int row = i;
            Future<?> future = executor.submit(() -> {
                data[row] = data[row].mod(p);
            });
            futures.add(future);
        }

        // 全ての行に対するタスクが完了するまで待機
        for (Future<?> future : futures) {
            try {
                future.get();
            } catch (InterruptedException | ExecutionException e) {
                e.printStackTrace();
            }
        }

        // スレッドプールをシャットダウン
        executor.shutdown();
	}


	public BigInteger[] getResult() {
		// for (int i = 0; i < MATRIX_SIZE; ++i) {
		// 	// for (int k = 0; k < MATRIX_SIZE; ++k) {
		// 	///	data[i] =  data[i].multiply(tmp[i][k]).mod(p);
		// 		data[i] =  data[i].mod(p);
		// 	//}
		// }
		c = data;
		return c;
	}

	private static BigInteger modPow(BigInteger base, BigInteger exponent, BigInteger mod) {
		if (base.equals(BigInteger.ZERO))
			return BigInteger.ZERO;
		else
			return base.modPow(exponent, mod);
	}

}

class PararllelSE2 {
	BigInteger[][] a;
	BigInteger[][] b;
	BigInteger[][] c;
	BigInteger [][] data;
	BigInteger p;
	int MATRIX_SIZE;
	int MINIMUM_THRESHOLD;

	private static final int
			 POOL_SIZE = Runtime.getRuntime().availableProcessors();

	private final ExecutorService exec = Executors.newFixedThreadPool(POOL_SIZE);

	PararllelSE2(BigInteger[][] a, BigInteger[][] b, BigInteger p, int d, int MINIMUM_THRESHOLD) {
		// assumption : a and b are both double[MATRIX_SIZE][MATRIX_SIZE]
		this.a = a;
		this.b = b;
		this.MATRIX_SIZE = d;
		this.p = p;
		this.data = new BigInteger[MATRIX_SIZE][MATRIX_SIZE];
		this.MINIMUM_THRESHOLD = MINIMUM_THRESHOLD;

		for (int i = 0; i < MATRIX_SIZE; ++i) {
			for (int j = 0; j < MATRIX_SIZE; ++j) {
				data[i][j] = BigInteger.ONE;
			}
		}
		System.out.println("POOL_SIZE:" + POOL_SIZE);
		//this.c = new double[MATRIX_SIZE][MATRIX_SIZE];
	}

	public void SERight() {
		//multiplyRecursive(0, 0, 0, 0, 0, 0, a.length);
		Future f = exec.submit(new SETaskRight(a, b, data, 0, 0, 0, 0, 0, 0, data.length));
		try {
			f.get();
			exec.shutdown();
		} catch (Exception e) {

		}
	}

	public void SELeft() {
		//multiplyRecursive(0, 0, 0, 0, 0, 0, a.length);
		Future f = exec.submit(new SETaskLeft(a, b, data, 0, 0, 0, 0, 0, 0, data.length));
		try {
			f.get();
			exec.shutdown();
		} catch (Exception e) {

		}
	}

	public BigInteger[][] getResult() {
		for (int i = 0; i < MATRIX_SIZE; ++i) {
			for (int j = 0; j < MATRIX_SIZE; ++j) {
				data[i][j] = data[i][j].mod(p);
			}
		}
		c = data;
		return c;
	}

	private static BigInteger modPow(BigInteger base, BigInteger exponent, BigInteger mod) {
		if (base.equals(BigInteger.ZERO))
			return BigInteger.ZERO;
		else
			return base.modPow(exponent, mod);
	}

	class SETaskRight implements Runnable{
		private BigInteger[][] a;
		private BigInteger[][] b;
		BigInteger [][] data;
		private int a_i, a_j, b_i, b_j, c_i, c_j, size;

		SETaskRight(BigInteger[][] a, BigInteger[][] b, BigInteger [][] data, int a_i, int a_j, int b_i, int b_j, int c_i, int c_j, int size) {
			this.a = a;
			this.b = b;
			this.data = data;
			this.a_i = a_i;
			this.a_j = a_j;
			this.b_i = b_i;
			this.b_j = b_j;
			this.c_i = c_i;
			this.c_j = c_j;
			this.size = size;
		}

		public void run() {
			//System.out.format("[%d,%d]x[%d,%d](%d)\n",a_i,a_j,b_i,b_j,size);
			int h = size/2;
			if (size <= MINIMUM_THRESHOLD) {
				for (int i = 0; i < size; ++i) {
					for (int j = 0; j < size; ++j) {
						for (int k = 0; k < size; ++k) {
							data[c_i+i][c_j+j] = data[c_i+i][c_j+j].multiply(modPow(a[a_i+i][a_j+k], b[b_i+k][b_j+j], p));

						}
						data[c_i+i][c_j+j].mod(p);
					}
				}
			} else {
				SETaskRight[] tasks = {
					new SETaskRight(a, b, data, a_i, a_j, b_i, b_j, c_i, c_j, h),
					new SETaskRight(a, b, data, a_i, a_j+h, b_i+h, b_j, c_i, c_j, h),

					new SETaskRight(a, b, data, a_i, a_j, b_i, b_j+h, c_i, c_j+h, h),
					new SETaskRight(a, b, data, a_i, a_j+h, b_i+h, b_j+h, c_i, c_j+h, h),

					new SETaskRight(a, b, data, a_i+h, a_j, b_i, b_j, c_i+h, c_j, h),
					new SETaskRight(a, b, data, a_i+h, a_j+h, b_i+h, b_j, c_i+h, c_j, h),

					new SETaskRight(a, b, data, a_i+h, a_j, b_i, b_j+h, c_i+h, c_j+h, h),
					new SETaskRight(a, b, data, a_i+h, a_j+h, b_i+h, b_j+h, c_i+h, c_j+h, h)
				};

				FutureTask[] fs = new FutureTask[tasks.length/2];

				for (int i = 0; i < tasks.length; i+=2) {
					fs[i/2] = new FutureTask(new SequentializerRight(tasks[i], tasks[i+1]), null);
					exec.execute(fs[i/2]);
				}
				for (int i = 0; i < fs.length; ++i) {
					fs[i].run();
				}
				try {
					for (int i = 0; i < fs.length; ++i) {
						fs[i].get();
					}
				} catch (Exception e) {

				}
			}
		}
	}

	class SETaskLeft implements Runnable{
		private BigInteger [][] a;
		private BigInteger [][] b;
		BigInteger [][] data;
		private int a_i, a_j, b_i, b_j, c_i, c_j, size;

		SETaskLeft(BigInteger [][] a, BigInteger [][] b, BigInteger [][] data, int a_i, int a_j, int b_i, int b_j, int c_i, int c_j, int size) {
			this.a = a;
			this.b = b;
			this.data = data;
			this.a_i = a_i;
			this.a_j = a_j;
			this.b_i = b_i;
			this.b_j = b_j;
			this.c_i = c_i;
			this.c_j = c_j;
			this.size = size;
		}

		public void run() {
			//System.out.format("[%d,%d]x[%d,%d](%d)\n",a_i,a_j,b_i,b_j,size);
			int h = size/2;
			if (size <= MINIMUM_THRESHOLD) {

				for (int i = 0; i < size; ++i) {
					for (int j = 0; j < size; ++j) {
						for (int k = 0; k < size; ++k) {
							data[c_i+i][c_j+j] = data[c_i+i][c_j+j].multiply(modPow(a[a_i + k][a_j + j], b[b_i + i][b_j + k], p));

						}
						data[c_i+i][c_j+j].mod(p);
					}
				}
			} else {
				SETaskLeft[] tasks = {
					new SETaskLeft(a, b, data, a_i, a_j, b_i, b_j, c_i, c_j, h),
					new SETaskLeft(a, b, data, a_i+h, a_j, b_i, b_j+h, c_i, c_j, h),

					new SETaskLeft(a, b, data, a_i, a_j, b_i+h, b_j, c_i+h, c_j, h),
					new SETaskLeft(a, b, data, a_i+h, a_j, b_i+h, b_j+h, c_i+h, c_j, h),

					new SETaskLeft(a, b, data, a_i, a_j+h, b_i, b_j, c_i, c_j+h, h),
					new SETaskLeft(a, b, data, a_i+h, a_j+h, b_i, b_j+h, c_i, c_j+h, h),

					new SETaskLeft(a, b, data, a_i, a_j+h, b_i+h, b_j, c_i+h, c_j+h, h),
					new SETaskLeft(a, b, data, a_i+h, a_j+h, b_i+h, b_j+h, c_i+h, c_j+h, h)
				};

				FutureTask[] fs = new FutureTask[tasks.length/2];

				for (int i = 0; i < tasks.length; i+=2) {
					fs[i/2] = new FutureTask(new SequentializerLeft(tasks[i], tasks[i+1]), null);
					exec.execute(fs[i/2]);
				}
				for (int i = 0; i < fs.length; ++i) {
					fs[i].run();
				}
				try {
					for (int i = 0; i < fs.length; ++i) {
						fs[i].get();
					}
				} catch (Exception e) {

				}
			}
		}
	}

	class SequentializerRight implements Runnable{
		private SETaskRight first, second;
		SequentializerRight(SETaskRight first, SETaskRight second) {
			this.first = first;
			this.second = second;
		}
		public void run() {
			first.run();
			second.run();
		}

	}

	class SequentializerLeft implements Runnable{
		private SETaskLeft first, second;
		SequentializerLeft(SETaskLeft first, SETaskLeft second) {
			this.first = first;
			this.second = second;
		}
		public void run() {
			first.run();
			second.run();
		}

	}



	public static BigInteger randomBInt(Random random, BigInteger modulus) {
		if (modulus.compareTo(new BigInteger("2")) < 0)
			throw new IllegalArgumentException("mod must be gleater or equalsto 2.");

		BigInteger res;
		BigInteger MAX = new BigInteger(Integer.toString(Integer.MAX_VALUE));

		if (modulus.compareTo(MAX) > 0) {
			BigInteger m = modulus;
			Vector<Integer> v = new Vector<Integer>();

			// Split biginteger into size of integer.
			while (!m.equals(BigInteger.ZERO)) {
				v.add(random.nextInt(m.mod(MAX).intValue()));
				m = m.divide(MAX);
			}

			res = BigInteger.ZERO;
			int cnt = 0;

			// Assembling the cut value again
			for (Integer b : v) {
				res = res.add(new BigInteger(b.toString()).multiply(MAX.pow(cnt)));
				cnt++;
			}

			return res;
		} else if (modulus.compareTo(BigInteger.ZERO) > 0) {
			// In case of mod < IntMAX
			do
				res = new BigInteger(Integer.toString(random.nextInt(modulus.intValue())));
			while (res.equals(BigInteger.ZERO));
			return res;
		} else
			throw new IllegalArgumentException("'mod' must be greater than 0.");
	}



}
