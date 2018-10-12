package com.cg.project.test;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import com.cg.project.bean.MathServices;
import com.cg.project.bean.MathServicesImpl;
import com.cg.project.exception.InvalidNuumberRangeException;
public class MathServicesTest {
	private static MathServices mathServices;
	@BeforeClass
	public static void setUPTestEnv() {
		mathServices = new MathServicesImpl();
	}
	@Before
	public void setUpMockDataForTest() {
	}
	
	@Test(expected=InvalidNuumberRangeException.class)
	public void testAddNumbersForFirstNOInvalid() throws InvalidNuumberRangeException {
		mathServices.addNums(-100, 200);
	}
	
	@Test(expected=InvalidNuumberRangeException.class)
	public void testAddNumbersForSecondNOInvalid() throws InvalidNuumberRangeException {
		mathServices.addNums(100, -200);
	}
	
	@Test()
	public void testAddNumbesrForBothValidNo() throws InvalidNuumberRangeException {
		int expectedAns=300;
		int actualAns=	mathServices.addNums(100, 200);
		Assert.assertEquals(expectedAns, actualAns);
	}
	@Test(expected=InvalidNuumberRangeException.class)
	public void testSubNumbersForFirstNOInvalid() throws InvalidNuumberRangeException {
		mathServices.subNums(-100, 200);
	}
	
	@Test(expected=InvalidNuumberRangeException.class)
	public void testSubNumbersForSecondNOInvalid() throws InvalidNuumberRangeException {
		mathServices.subNums(100, -200);
	}
	
	@Test()
	public void testSubNumbesrForBothValidNo() throws InvalidNuumberRangeException {
		int expectedAns=100;
		int actualAns=	mathServices.subNums(200, 100);
		Assert.assertEquals(expectedAns, actualAns);
	}
	@Test(expected=InvalidNuumberRangeException.class)
	public void testMulNumbersForFirstNOInvalid() throws InvalidNuumberRangeException {
		mathServices.multiNums(-100, 200);
	}
	
	@Test(expected=InvalidNuumberRangeException.class)
	public void testMulNumbersForSecondNOInvalid() throws InvalidNuumberRangeException {
		mathServices.multiNums(100, -200);
	}
	
	@Test()
	public void testMulNumbesrForBothValidNo() throws InvalidNuumberRangeException {
		int expectedAns=200;
		int actualAns=	mathServices.multiNums(10, 20);
		Assert.assertEquals(expectedAns, actualAns);
	}
	
	@After
	public void tearDownMockDataForTest() {
		System.out.println("tearDownMockDataForTest()");
	}
	
	
	@AfterClass
	public static void tearDownTestEnv() {
		mathServices=null;
	}
}
