Monolithic implementation of the Range Minimum Query (RMQ) with _Cartesian trees_ based on [this](https://github.com/birc-stormtroopers/rmq).

Since a _Sparse table_ is needed, it also contains its implementation and uses it for testing.

**Note** that both implementations return the _value_ of the minimum in the range, and not the _index_ of the minimum.

It is only intended as a reference. Tests can be run with `cargo test`.