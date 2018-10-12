package com.cg.project.bean;
import com.cg.project.exception.InvalidNuumberRangeException;
public class MathServicesImpl implements MathServices{
	@Override
	public int addNums(int n1, int n2) throws InvalidNuumberRangeException{
		if(n1<0||n2<0)
			throw new InvalidNuumberRangeException("Enter a number greater than zero");
		return n1+n2;
	}

	@Override
	public int subNums(int n1, int n2) throws InvalidNuumberRangeException{
		if(n1<0||n2<0)
			throw new InvalidNuumberRangeException("Enter a number greater than zero");
		return n1-n2;
	}

	@Override
	public int multiNums(int n1, int n2) throws InvalidNuumberRangeException{
		if(n1<0||n2<0)
			throw new InvalidNuumberRangeException("Enter a number greater than zero");
		return n1*n2;
	}
	

}
