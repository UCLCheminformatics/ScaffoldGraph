## Contributing

Contributions to ScaffoldGraph will most likely fall into the following categories:

1. Implementing a new Feature:
    * New Features that fit into the scope of this package will be accepted. If you are unsure about the 
      idea/design/implementation, feel free to post an issue.
2. Fixing a Bug:
    * Bug fixes are welcomed, please send a Pull Request each time a bug is encountered. When sending a Pull
      Request please provide a clear description of the encountered bug. If unsure feel free to post an issue

Please send Pull Requests to: 
http://github.com/UCLCheminformatics/scaffoldgraph

### Testing

ScaffoldGraphs testing is located under `test/`. Run all tests using:

```
$ python setup.py test
```

or run an individual test: `pytest --no-cov tests/core`

When contributing new features please include appropriate test files

### Continuous Integration

ScaffoldGraph uses Travis CI for continuous integration
