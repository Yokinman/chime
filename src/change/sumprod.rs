//! ... https://www.desmos.com/calculator/ko5owr56jx

#[cfg(test)]
mod tests {
	use crate as chime;
	use crate::Flux;
	use crate::poly::{Poly, Roots};

	#[derive(Flux)]
	struct Test {
		#[basis]
		value: f64,
		#[change(add_per(chime::time::SEC))]
		add: f64,
		#[change(mul_per(chime::time::SEC))]
		mul: f64,
	}
	
	#[derive(Flux)]
	struct Test2 {
		#[basis]
		value: f64,
		#[change(add_per(chime::time::SEC))]
		add: Test,
	}
	
	#[test]
	fn sumprod() {
		let a = Test { value: 1., add: 4., mul: 0.5 }.to_poly();
		assert_eq!(a.eval(0.), 1.);
		assert_eq!(a.eval(1.), 1.9426950408889634);
		assert_eq!(a.eval(2.), 2.414042561333445);
		assert_eq!(a.eval(3.), 2.649716321555686);
		assert_eq!(a.eval(f64::INFINITY), 2.8853900817779268);
		assert_eq!(a.eval(f64::NEG_INFINITY), f64::NEG_INFINITY);
		let a = Test { value: 1., add: 4., mul: 2. }.to_poly();
		assert_eq!(a.eval(0.), 1.);
		assert_eq!(a.eval(1.), 13.541560327111707);
		assert_eq!(a.eval(2.), 38.62468098133512);
		assert_eq!(a.eval(3.), 88.79092228978195);
		assert_eq!(a.eval(f64::INFINITY), f64::INFINITY);
		assert_eq!(a.eval(f64::NEG_INFINITY), -11.541560327111707);
		
		for root in (a - symb_poly::Invar(chime::constant::Constant(5.))).roots() {
			let val = a.eval(root);
			assert!((val - 5.).abs() < 1e-12, "{:?}", val);
		}
		
		let b = Test2 {
			value: 1.,
			add: Test { value: 1., add: 4., mul: 2. }
		}.to_poly();
		assert_eq!(b.eval(0.), 1.);
		assert_eq!(b.eval(1.), 7.552086561822119);
		assert_eq!(b.eval(2.), 32.19782001257806);
		assert_eq!(b.eval(3.), 93.03084724120166);
		assert_eq!(b.eval(-1.), 3.494736882644794);
		assert_eq!(b.eval(-2.), 10.512885487523043);
		
		for root in (b - symb_poly::Invar(chime::constant::Constant(5.))).roots() {
			let val = b.eval(root);
			assert!((val - 5.).abs() < 1e-12, "{:?}", (val, root));
		}
		
		let c = Test { value: 1., add: 4., mul: 1. }.to_poly();
		assert_eq!(c.eval(0.), 1.);
		assert_eq!(c.eval(1.), 5.);
		assert_eq!(c.eval(2.), 9.);
		assert_eq!(c.eval(3.), 13.);
		
		// for root in (c - symb_poly::Invar(chime::constant::Constant(5.))).roots() {
		// 	let val = c.eval(root);
		// 	assert!((val - 5.).abs() < 1e-12, "{:?}", val);
		// }
		
		let d = Test2 {
			value: 1.,
			add: Test { value: 1., add: 4., mul: 1. }
		}.to_poly();
		assert_eq!(d.eval(0.), 1.);
		assert_eq!(d.eval(1.), 4.);
		assert_eq!(d.eval(2.), 11.);
		assert_eq!(d.eval(3.), 22.);
		
		todo!("re-add Roots exception for mul=1 case");
		// for root in (d - symb_poly::Invar(chime::constant::Constant(5.))).roots() {
		// 	let val = d.eval(root);
		// 	assert!((val - 5.).abs() < 1e-12, "{:?}", val);
		// }
	}
}