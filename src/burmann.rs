
fn partial_partition_array( partition :&mut Vec<u64>,  vec_partition: &mut Vec<Vec<u64>>, max_i : u64, k :u64){
    // for i in 1..max_i{
    //     partition[i] = 0;
    // }
    // parition[0] = k-max_i;
    // vec_partition.push(partition);
    for i in 2..max_i{
        let index = (i-1) as usize;
        while partition[0] >= i{
            partition[0] -=i;
            let prev= partition[0];
            partition[index] +=1;
            vec_partition.push(partition.to_vec());
            partial_partition_array(partition, vec_partition, i, partition[0]);
            partition[0] = prev;
        }
        partition[index] = 0;
        partition[0] = k;
    }
    //resetting
    for i in 0..max_i as usize{
        partition[i] = 0;
    }
}

pub fn integer_partition_array(k:u64) -> Vec<Vec<u64>>{
    let mut partition = vec![0_u64; k as usize];
    let mut vec_partition= Vec::<Vec::<u64>>::new();
    let last_index = (k-1) as usize;
    partition[last_index] = 1;
    vec_partition.push(partition.to_vec());
    partition[last_index] = 0;
    partition[0] = k;
    vec_partition.push(partition.to_vec());
    partial_partition_array(&mut partition, &mut vec_partition, k,k);
    vec_partition
}