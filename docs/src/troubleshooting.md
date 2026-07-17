# Troubleshooting

## Sandia network SSL/TLS errors

If you are on Sandia's network and run into SSL/TLS certificate errors when
installing or fetching packages, see `README-sandia.md` in the repository for a
complete setup guide.

## Neural-network model aborts asking for PyCall

The neural-network Operator-Inference model is an optional package extension
that needs Python (PyTorch via PyCall). If you request it without PyCall
installed, Norma aborts with a message telling you to install it. See
[Installation](installation.md).

## A required input key is missing

Most input keys are required and are read directly; a missing one aborts with a
message naming the key. Check the section for that key in the
[Input File Reference](reference/index.md), where each key is marked required or
optional with its default.

## License

Norma is licensed under the BSD 3-Clause License. See `LICENSE.md` in the
repository for details.
