#include "BuildingBlocks/aux-protocols.h"
#include "utils/emp-tool.h"
#include <iostream>
using namespace sci;
using namespace std;

int party, port = 8000, dim = 1 << 16;
string address = "127.0.0.1";
NetIO *io;
OTPack<NetIO> *otpack;
AuxProtocols *aux;


// test_ring_aux_protocols
int main(int argc, char **argv)
{
    cout<<"开始了uu们"<<endl;
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("d", dim, "Size of vector");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.parse(argc, argv);

    io = new NetIO(party == 1 ? nullptr : address.c_str(), port);
    otpack = new OTPack<NetIO>(io, party);

    aux = new AuxProtocols(party, io, otpack);

    // -->test_wrap_computation
    int bw_x = 32;
    PRG128 prg;
    uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));

    uint64_t *x = new uint64_t[dim];
    uint8_t *y = new uint8_t[dim];

    prg.random_data(x, dim * sizeof(uint64_t)); // 生成随机数
    for (int i = 0; i < dim; i++)
    {
        x[i] = x[i] & mask_x;
    }// 感觉是把最高位当成符号,取其他位了

    // -->-->wrap_computation
    assert(bw_x <= 64);
    uint64_t mask = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));

    uint64_t *tmp_x = new uint64_t[dim];

    for (int i = 0; i < dim; i++) {
    if (party == sci::ALICE)
        tmp_x[i] = x[i] & mask;
    else
        tmp_x[i] = (mask - x[i]) & mask; // 2^{bw_x} - 1 - x[i]
    }

    aux->mill->compare(y, tmp_x, dim, bw_x, true); // computing greater_than


    if (party == ALICE) {
        io->send_data(x, dim * sizeof(uint64_t));
        io->send_data(y, dim * sizeof(uint8_t));
    } else {
        uint64_t *x0 = new uint64_t[dim];
        uint8_t *y0 = new uint8_t[dim];
        io->recv_data(x0, dim * sizeof(uint64_t));
        io->recv_data(y0, dim * sizeof(uint8_t));

        for (int i = 0; i < dim; i++) {
        assert((x0[i] > (mask_x - x[i])) == (y0[i] ^ y[i]));
        }
        cout << "Wrap Computation Tests passed" << endl;

        delete[] x0;
        delete[] y0;
    }
    return 0;
}
/*
整个过程
*/