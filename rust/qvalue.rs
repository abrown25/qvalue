use std::io::BufferedReader;
use std::io::File;

fn main(){

    let path = Path::new("vQTLresults.txt");
    let mut file = BufferedReader::new(File::open(&path));
    let _ = file.read_line();

    let mut p_val = Vec::new();
    for line in file.lines()
    {
        let v = from_str::<f32>(line.unwrap().as_slice().words().nth(3).unwrap().as_slice());
        match v 
        {
            Some(v) => {
                if v >= 0.0 && v <= 1.0 { p_val.push(v); }
            },
            None => ()
        }
    }
    
    let mut order_index: Vec<uint> = range(0, p_val.len()).collect();
    order_index.sort_by(|&a, &b| if p_val.get(a) > p_val.get(b) { Less } else if p_val.get(a)==p_val.get(b) { Equal } else { Greater });

    let lambda = Vec::from_fn(19, |x| 0.9 - (x as f32)*0.05);
    let mut pi0 :Vec<f32> = Vec::with_capacity(lambda.len());

    let mut val_iter = lambda.iter();
    let mut value = val_iter.next().unwrap();
    let mut number = 0.0;
    
    let mut ord_iter = order_index.iter();
    for place in ord_iter
    {
        if p_val.get(*place) < value {
            pi0.push(number / (1.0 - *value) / (p_val.len() as f32));
            let v = val_iter.next();
            match v {
                Some(v) => { 
                    value = v;
                }
                None => {
                    break;
                } 
            }
        }
        number += 1.0;
    }

    if ord_iter.next()==None{
        pi0.push(1.0);
    }

    println!("{} {}", pi0, lambda);
}
