package com.cg.project.client;

import com.cg.project.bean.MathServices;
import com.cg.project.bean.MathServicesImpl;
import com.cg.project.exception.InvalidNuumberRangeException;

public class MainClass {

	public static void main(String[] args) {
		try{
			MathServices services = new MathServicesImpl();
			System.out.println("Addition "+services.addNums(5, 20));
			System.out.println("Substraction "+services.subNums(3, 9));
			System.out.println("Multiplication "+services.multiNums(5, -5));
		}
		catch(InvalidNuumberRangeException e){
			System.out.println("Number should not be less than zero");
			e.printStackTrace();
		}

	}

}
