const {
	acos,
	acosh,
	asinh,
	cos,
	cosh,
	sinh,
	sqrt,
	abs,
	sign,
	PI
} = Math;

/**
 * Находит рациональные корни уравнения x^3+ax^2+bx+c=0.
 */
function cube_viet(a, b, c){
	const Q = (a**2 - 3*b)/9;
	const R = (2*a**3 - 9*a*b + 27*c)/54;
	const Q3 = Q**3;
	const Q15 = sqrt(abs(Q3));
	const Q05 = sqrt(abs(Q));
	const S = Q3 - R**2;
	//console.log(S,Q);
	if(S>0){
		const varphi = R/Q15;
		const phi = acos(varphi)/3;
		const step = 2*PI/3;
		
		let x = [-step, 0, step].map((step)=>(
			-2*Q05*cos(phi+step) - a/3
		));
		
		return x.sort();
	}
	else if(S<0){
		if(Q!==0){
			const varphi = abs(R)/Q15;
			if(Q>0){
				const phi = acosh(varphi)/3;
				
				let x = -2 * sign(R) * Q05 * cosh(phi) - a/3;
				return [x];
			}
			else{
				const phi = asinh(varphi)/3;
				let x = -2 * Q05 * sinh(phi) - a/3;
				return [x];
			}
		}
		else{
			console.log(c,a);
			let x = -(Math.cbrt(c - (a**3)/27)) - a/3;
			return [x];
		}
	}
	else if(S===0){
		const R_3 = R**(1/3);
		let x = [
			-2*R_3-a/3,
			R_3-a/3
		]
		return x.sort();
	}
	return [];
}

/**
 * Находит рациональные корни уравнения ax^2+bx+c =0;
 */
function square(a, b, c){
	if(a){
		const D = b**2 - 4 *a*c;
		if(D>0){
			let step = sqrt(D);
			let x = [-step, step].map((step)=>(
				(-b + step)/2/a
			));
			return x;
		}
		else if(D===0){
			let x = -b/2/a;
			return [x];
		}
		else{
			return [];
		}
	}
	else{
		
		return linear(b, c);
	}
	return [];
}

/**
 * Находит конечные корни уравнения ax+b=0
 */
function linear(a, b){
	if(a){
		let x = -b/a;
		return [x];
	}
	return [];
}

/**
 * Находит рациональные корни уравнения ax^3+bx^2+cx+t=0;
 */
function cube(a, b, c, d){
	a=a||0; b=b||0; c=c||0; d=d||0;
	if(a){
		return cube_viet(b/a, c/a, d/a);
	}
	else{
		return square(b, c, d);
	}
}

module.exports = cube;