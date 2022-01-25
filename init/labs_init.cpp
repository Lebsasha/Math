// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com


//labs_init::labs_init()
//{
//    for (int i = 2; i <= 7; ++i)
//    {
//        path path = ".";
//        path /= std::to_string(i) +"_Lab";
//        assert(exists(path) && is_directory(path));
//        path/="output";
//        if (!exists(path))
//            create_directory(path);
//    }
//
//    path lab_path = "./1_Lab/input";
//    if (!exists(lab_path))
//    {
//        create_directory(lab_path);
//        double p[] = {1, 2, 1, -1, -2, 2, 0, 1, 1};
//        Matrix_SLE A (p, 3, 3);
//        A.write_to_file((lab_path / "data.matr_l").c_str());
//    }
//}