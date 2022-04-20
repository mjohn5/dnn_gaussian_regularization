# dnn_gaussian_regularization

contains codes (developed by Sujit Vettam) used for Gaussian regularization presented in John, Vettam and Wu (2022): arXiv: 2204.03123

https://arxiv.org/pdf/2204.03123.pdf

1) codes for DNN architecture applied on MNIST dataset is included (the same DNN architecture was applied on FMNIST and RCV1 datasets)

2) codes for CNN architecture applied on SVHN and CIFAR-100 datasets are also included (the same CNN architecture was applied on CIFAR-10 and ImagNet datasets mentioned in the paper). CIFAR-10 codes can be obtained by changing 'cifar100' to 'cifar10' in the line '(X_train, y_train), (X_test, y_test) = tf.keras.datasets.cifar100.load_data()' within the CIFAR-100 codes.
